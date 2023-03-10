/* This file is part of the KDE project
 * SPDX-FileCopyrightText: 2016 Julian Thijssen <julianthijssen@gmail.com>
 * SPDX-FileCopyrightText: 2023 L. E. Segovia <amy@amyspark.me>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "kis_opengl_shader_loader.h"

#include "opengl/kis_opengl.h"
#include "kis_config.h"

#include <config-ocio.h>

#include <QFile>
#include <QMessageBox>
#include <KLocalizedString>

#define PROGRAM_VERTEX_ATTRIBUTE 0
#define PROGRAM_TEXCOORD_ATTRIBUTE 1

// Mapping of uniforms to uniform names
std::map<Uniform, const char *> KisShaderProgram::names = {
   {ModelViewProjection, "modelViewProjection"},
   {TextureMatrix, "textureMatrix"},
   {ViewportScale, "viewportScale"},
   {TexelSize, "texelSize"},
   {Texture0, "texture0"},
   {Texture1, "texture1"},
   {FixedLodLevel, "fixedLodLevel"},
   {FragmentColor, "fragColor"}
};

/**
 * Generic shader loading function that will compile a shader program given
 * a vertex shader and fragment shader resource path. Extra code can be prepended
 * to each shader respectively using the header parameters.
 *
 * @param vertPath Resource path to a vertex shader
 * @param fragPath Resource path to a fragment shader
 * @param vertHeader Extra code which will be prepended to the vertex shader
 * @param fragHeader Extra code which will be prepended to the fragment shader
 */
KisShaderProgram *KisOpenGLShaderLoader::loadShader(QString vertPath, QString fragPath,
                                                    QByteArray vertHeader, QByteArray fragHeader)
{
    bool result;

    KisShaderProgram *shader = new KisShaderProgram();

    // Load vertex shader
    QByteArray vertSource;

    if (KisOpenGL::hasOpenGLES()) {
        vertSource.append("#version 300 es\n");
    } else {
#ifdef Q_OS_MACOS
        vertSource.append(KisOpenGL::hasOpenGL3() ? "#version 150 core\n" : "#version 120\n");
        vertSource.append("#define texture2D texture\n");
        vertSource.append("#define texture3D texture\n");
#else
        vertSource.append(KisOpenGL::supportsLoD() ? "#version 130\n" : "#version 120\n");
#endif
    }
    vertSource.append(vertHeader);
    QFile vertexShaderFile(":/" + vertPath);
    vertexShaderFile.open(QIODevice::ReadOnly);
    vertSource.append(vertexShaderFile.readAll());

    result = shader->addShaderFromSourceCode(QOpenGLShader::Vertex, vertSource);
    if (!result)
        throw ShaderLoaderException(QString("%1: %2 - Cause: %3").arg("Failed to add vertex shader source from file", vertPath, shader->log()));

    // Load fragment shader
    QByteArray fragSource;

    if (KisOpenGL::hasOpenGLES()) {
        fragSource.append("#version 300 es\n");
        if (KisOpenGL::supportsLoD()) {
            fragSource.append("#extension GL_EXT_shader_texture_lod : enable\n");
        }
        fragSource.append(
            "precision mediump float;\n"
            "precision mediump sampler3D;\n");

        // OpenColorIO doesn't support OpenGL ES.
        fragSource.append("#define texture2D texture\n");
        fragSource.append("#define texture3D texture\n");
        if (KisOpenGL::supportsLoD()) {
            fragSource.append(
                "#if __VERSION__ < 300\n"
                "#define textureLod texture2DLodEXT\n"
                "#endif\n"
            );
        }
    } else {
#ifdef Q_OS_MACOS
        fragSource.append(KisOpenGL::hasOpenGL3() ? "#version 150 core\n" : "#version 120\n");
        fragSource.append("#define texture2D texture\n");
        fragSource.append("#define texture3D texture\n");
#else
        fragSource.append(KisOpenGL::supportsLoD() ? "#version 130\n" : "#version 120\n");
#endif
    }
    fragSource.append(fragHeader);
    QFile fragmentShaderFile(":/" + fragPath);
    fragmentShaderFile.open(QIODevice::ReadOnly);
    fragSource.append(fragmentShaderFile.readAll());

    result = shader->addShaderFromSourceCode(QOpenGLShader::Fragment, fragSource);
    if (!result)
        throw ShaderLoaderException(QString("%1: %2 - Cause: %3").arg("Failed to add fragment shader source from file", fragPath, shader->log()));

    // Bind attributes
    shader->bindAttributeLocation("a_vertexPosition", PROGRAM_VERTEX_ATTRIBUTE);
    shader->bindAttributeLocation("a_textureCoordinate", PROGRAM_TEXCOORD_ATTRIBUTE);

    // Link
    result = shader->link();
    if (!result)
        throw ShaderLoaderException(QString("Failed to link shader: ").append(vertPath));

    Q_ASSERT(shader->isLinked());

    return shader;
}

/**
 * Specific display shader loading function. It adds the appropriate extra code
 * to the fragment shader depending on what is available on the target machine.
 * Additionally, it picks the appropriate shader files depending on the availability
 * of OpenGL3.
 */
KisShaderProgram *KisOpenGLShaderLoader::loadDisplayShader(QSharedPointer<KisDisplayFilter> displayFilter, bool useHiQualityFiltering)
{
    QByteArray fragHeader;

    if (KisOpenGL::supportsLoD()) {
        fragHeader.append("#define DIRECT_LOD_FETCH\n");
        if (useHiQualityFiltering) {
            fragHeader.append("#define HIGHQ_SCALING\n");
        }
    }

    // If we have an OCIO display filter and it contains a function we add
    // it to our shader header which will sit on top of the fragment code.
    bool haveDisplayFilter = displayFilter && !displayFilter->program().isEmpty();
    if (haveDisplayFilter) {
        fragHeader.append("#define USE_OCIO\n");
#ifdef HAVE_OCIO_V2
        fragHeader.append("#define USE_OCIO_V2\n");
#endif
        fragHeader.append(displayFilter->program().toLatin1());
    }

    QString vertPath, fragPath;
    // Select appropriate shader files
    if (KisOpenGL::supportsLoD()) {
        vertPath = "matrix_transform.vert";
        fragPath = "highq_downscale.frag";
    } else {
        vertPath = "matrix_transform_legacy.vert";
        fragPath = "simple_texture_legacy.frag";
    }

    KisShaderProgram *shader = loadShader(vertPath, fragPath, QByteArray(), fragHeader);

    return shader;
}

/**
 * Specific checker shader loading function. It picks the appropriate shader
 * files depending on the availability of OpenGL3 on the target machine.
 */
KisShaderProgram *KisOpenGLShaderLoader::loadCheckerShader()
{
    QString vertPath, fragPath;
    // Select appropriate shader files
    if (KisOpenGL::supportsLoD()) {
        vertPath = "matrix_transform.vert";
        fragPath = "simple_texture.frag";
    } else {
        vertPath = "matrix_transform_legacy.vert";
        fragPath = "simple_texture_legacy.frag";
    }

    KisShaderProgram *shader = loadShader(vertPath, fragPath, QByteArray(), QByteArray());

    return shader;
}

/**
 * Specific uniform shader loading function. It picks the appropriate shader
 * files depending on the availability of OpenGL3 on the target machine.
 */
KisShaderProgram *KisOpenGLShaderLoader::loadSolidColorShader()
{
    QString vertPath, fragPath;
    // Select appropriate shader files
    if (KisOpenGL::supportsLoD()) {
        vertPath = "solid_color.vert";
        fragPath = "solid_color.frag";
    } else {
        vertPath = "solid_color_legacy.vert";
        fragPath = "solid_color_legacy.frag";
    }

    KisShaderProgram *shader = loadShader(vertPath, fragPath, QByteArray(), QByteArray());

    return shader;
}
