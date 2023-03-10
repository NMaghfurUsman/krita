/*
 *  SPDX-FileCopyrightText: 2016 Dmitry Kazakov <dimula73@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "VideoExportOptionsDialog.h"
#include "ui_video_export_options_dialog.h"

#include <KoID.h>

#include <ksharedconfig.h>
#include <kconfiggroup.h>
#include "VideoHDRMetadataOptionsDialog.h"
#include "KisHDRMetadataOptions.h"


struct KisVideoExportOptionsDialog::Private
{
    Private(ContainerType _containerType)
        : containerType(_containerType)
    {
        switch (containerType) {
            case ContainerType::OGV:
                codecs << KoID("libtheora", i18nc("theora codec name, check simplescreenrecorder for standard translations", "Theora"));
                break;
            case ContainerType::WEBM:
                codecs << KoID("libvpx", i18nc("VP9 codec name", "VP9"));
                break;
            case ContainerType::GIF:
                codecs << KoID("gif", i18nc("GIF codec name", "GIF"));
                break;
            case ContainerType::APNG:
                codecs << KoID("apng", i18nc("APNG codec name", "APNG"));
                break;
            case ContainerType::WEBP:
                codecs << KoID("libwebp", i18nc("WEBP codec name", "WEBP"));
                break;
            default:
                codecs << KoID("libx264", i18nc("h264 codec name, check simplescreenrecorder for standard translations", "H.264, MPEG-4 Part 10"));
                codecs << KoID("libx265", i18nc("h265 codec name, check simplescreenrecorder for standard translations", "H.265, MPEG-H Part 2 (HEVC)"));
                break;
        }

        presets << KoID("ultrafast", i18nc("h264 preset name, check simplescreenrecorder for standard translations", "ultrafast"));
        presets << KoID("superfast", i18nc("h264 preset name, check simplescreenrecorder for standard translations", "superfast"));
        presets << KoID("veryfast", i18nc("h264 preset name, check simplescreenrecorder for standard translations", "veryfast"));
        presets << KoID("faster", i18nc("h264 preset name, check simplescreenrecorder for standard translations", "faster"));
        presets << KoID("fast", i18nc("h264 preset name, check simplescreenrecorder for standard translations", "fast"));
        presets << KoID("medium", i18nc("h264 preset name, check simplescreenrecorder for standard translations", "medium"));
        presets << KoID("slow", i18nc("h264 preset name, check simplescreenrecorder for standard translations", "slow"));
        presets << KoID("slower", i18nc("h264 preset name, check simplescreenrecorder for standard translations", "slower"));
        presets << KoID("veryslow", i18nc("h264 preset name, check simplescreenrecorder for standard translations", "veryslow"));
        presets << KoID("placebo", i18nc("h264 preset name, check simplescreenrecorder for standard translations", "placebo"));

        profilesH264 << KoID("baseline", i18nc("h264 profile name, check simplescreenrecorder for standard translations", "baseline"));
        profilesH264 << KoID("main", i18nc("h264 profile name, check simplescreenrecorder for standard translations", "main"));
        profilesH264 << KoID("high", i18nc("h264 profile name, check simplescreenrecorder for standard translations", "high"));
        profilesH264 << KoID("high10", i18nc("h264 profile name, check simplescreenrecorder for standard translations", "high10"));
        profilesH264 << KoID("high422", i18nc("h264 profile name, check simplescreenrecorder for standard translations", "high422"));
        profilesH264 << KoID("high444", i18nc("h264 profile name, check simplescreenrecorder for standard translations", "high444"));

        profilesH265 << KoID("main", i18nc("h264 profile name, check simplescreenrecorder for standard translations", "main"));
        profilesH265 << KoID("main10", i18nc("h264 profile name, check simplescreenrecorder for standard translations", "main10 (HDR)"));

        // TODO: add "none" tune option
        tunesH264 << KoID("film", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "film"));
        tunesH264 << KoID("animation", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "animation"));
        tunesH264 << KoID("grain", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "grain"));
        tunesH264 << KoID("stillimage", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "stillimage"));
        tunesH264 << KoID("psnr", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "psnr"));
        tunesH264 << KoID("ssim", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "ssim"));
        tunesH264 << KoID("fastdecode", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "fastdecode"));
        tunesH264 << KoID("zerolatency", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "zerolatency"));

        tunesH265 << KoID("none", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "none"));
        tunesH265 << KoID("animation", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "animation"));
        tunesH265 << KoID("grain", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "grain"));
        tunesH265 << KoID("psnr", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "psnr"));
        tunesH265 << KoID("ssim", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "ssim"));
        tunesH265 << KoID("fastdecode", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "fastdecode"));
        tunesH265 << KoID("zero-latency", i18nc("h264 tune option name, check simplescreenrecorder for standard translations", "zero-latency"));
        
        predAPNG << KoID("none", i18nc("apng prediction option name", "none"));
        predAPNG << KoID("sub", i18nc("apng prediction option name", "sub"));
        predAPNG << KoID("up", i18nc("apng prediction option name", "up"));
        predAPNG << KoID("avg", i18nc("apng prediction option name", "avg"));
        predAPNG << KoID("paeth", i18nc("apng prediction option name", "paeth"));
        predAPNG << KoID("mixed", i18nc("apng prediction option name", "mixed"));

        presetWEBP << KoID("default", i18nc("webp preset option name", "default"));        
        presetWEBP << KoID("none", i18nc("webp preset option name", "none"));
        presetWEBP << KoID("drawing", i18nc("webp preset option name", "drawing"));        
        presetWEBP << KoID("icon", i18nc("webp preset option name", "icon"));
        presetWEBP << KoID("photo", i18nc("webp preset option name", "photo"));
        presetWEBP << KoID("picture", i18nc("webp preset option name", "picture"));
        presetWEBP << KoID("text", i18nc("webp preset option name", "text"));
        
        palettegenStatsMode << KoID("full", i18nc("paletegen status mode option name", "Global/Full"));
        palettegenStatsMode << KoID("diff", i18nc("paletegen status mode option name", "Difference"));
        palettegenStatsMode << KoID("single", i18nc("paletegen status mode option name", "Per Single Frame"));
    
        paletteuseDither << KoID("none", i18nc("paleteuse dither option name", "none"));
        paletteuseDither << KoID("bayer", i18nc("paleteuse dither option name", "bayer"));
        paletteuseDither << KoID("floyd_steinberg", i18nc("paleteuse dither option name", "floyd_steinberg"));
        paletteuseDither << KoID("heckbert", i18nc("paleteuse dither option name", "heckbert"));
        paletteuseDither << KoID("sierra2", i18nc("paleteuse dither option name", "sierra2"));
        paletteuseDither << KoID("sierra2_4a", i18nc("paleteuse dither option name", "sierra2_4a"));
        
        paletteuseDiffMode << KoID("none", i18nc("paleteuse diff mode option name", "none"));
        paletteuseDiffMode << KoID("rectangle", i18nc("paleteuse diff mode option name", "rectangle"));
    }

    QVector<KoID> codecs;
    QVector<KoID> presets;
    QVector<KoID> profilesH264;
    QVector<KoID> profilesH265;
    
    QVector<KoID> predAPNG;
    
    QVector<KoID> palettegenStatsMode;
    QVector<KoID> paletteuseDither;
    QVector<KoID> paletteuseDiffMode;

    QVector<KoID> tunesH264;
    QVector<KoID> tunesH265;
    QVector<KoID> presetWEBP;
    bool supportsHDR = false;
    ContainerType containerType;

    QString currentCustomLine;

    KisHDRMetadataOptions hdrMetadataOptions;
};

void populateComboWithKoIds(QComboBox *combo, const QVector<KoID> &ids, int defaultIndex)
{
    Q_FOREACH (const KoID &id, ids) {
        combo->insertItem(combo->count(), id.name());
    }
    combo->setCurrentIndex(defaultIndex);
}

KisVideoExportOptionsDialog::KisVideoExportOptionsDialog(ContainerType containerType, QWidget *parent)
    : KisConfigWidget(parent),
      ui(new Ui::VideoExportOptionsDialog),
      m_d(new Private(containerType))
{
    ui->setupUi(this);

    ui->intCRFH264->setRange(0, 51);
    ui->intCRFH264->setValue(28);

    ui->intCRFH265->setRange(0, 51);
    ui->intCRFH265->setValue(28);

    populateComboWithKoIds(ui->cmbPresetH264, m_d->presets, 5);
    populateComboWithKoIds(ui->cmbPresetH265, m_d->presets, 5);

    populateComboWithKoIds(ui->cmbProfileH264, m_d->profilesH264, 0);
    populateComboWithKoIds(ui->cmbProfileH265, m_d->profilesH265, 0);

    populateComboWithKoIds(ui->cmbTuneH264, m_d->tunesH264, 0);
    populateComboWithKoIds(ui->cmbTuneH265, m_d->tunesH265, 0);

    ui->intBitrate->setRange(10, 50000);
    ui->intBitrate->setValue(5000);
    ui->intBitrate->setSuffix(i18nc("kilo-bits-per-second, video bitrate suffix", "kbps"));

    ui->gifAutoPalette->setChecked(true);
    ui->gifReserveTransparent->setChecked(true);
    ui->gifLoop->setChecked(true);
    ui->gifTransDiff->setChecked(true);
    
    populateComboWithKoIds(ui->cmbPalettegenStatsModeGIF, m_d->palettegenStatsMode, 0);
    populateComboWithKoIds(ui->cmbPaletteuseDitherGIF, m_d->paletteuseDither, 5);
    populateComboWithKoIds(ui->cmbPaletteuseDiffModeGIF, m_d->paletteuseDiffMode, 0);
 
    ui->intPaletteuseBayerScaleGIF->setRange(0, 5);
    ui->intPaletteuseBayerScaleGIF->setValue(2);
    
    ui->apngLoop->setChecked(true);
    
    populateComboWithKoIds(ui->cmbPredAPNG, m_d->predAPNG, 0);
    
    ui->intCompressWEBP->setRange(0, 6);
    ui->intCompressWEBP->setValue(4);

    ui->intQscaleWEBP->setRange(0, 100);
    ui->intQscaleWEBP->setValue(75);
    
    populateComboWithKoIds(ui->cmbPresetWEBP, m_d->presetWEBP, 0);
    
    ui->webpLoop->setChecked(true);
    
    populateComboWithKoIds(ui->cmbCodec, m_d->codecs, 0);
    connect(ui->cmbCodec, SIGNAL(currentIndexChanged(int)), SLOT(slotCodecSelected(int)));
    slotCodecSelected(0);

    // TODO: temporarily hidden! Some combinations of 'tune' and
    //       'profile' options make ffmpeg generate empty file.
    //       We should not let the user shoot into his own foot!
    ui->cmbTuneH264->setVisible(false);
    ui->lblTuneH264->setVisible(false);

    ui->cmbTuneH265->setVisible(false);
    ui->lblTuneH265->setVisible(false);

    slotCustomLineToggled(false);
    connect(ui->chkCustomLine, SIGNAL(toggled(bool)), SLOT(slotCustomLineToggled(bool)));
    connect(ui->txtCustomLine, SIGNAL(editingFinished()), SLOT(slotSaveCustomLine()));
    connect(ui->btnResetCustomLine, SIGNAL(clicked()), SLOT(slotResetCustomLine()));

    connect(ui->chkUseHDRMetadata, SIGNAL(toggled(bool)),
            ui->btnHdrMetadata, SLOT(setEnabled(bool)));
    connect(ui->cmbProfileH265,
            SIGNAL(currentIndexChanged(int)),
            SLOT(slotH265ProfileChanged(int)));
    slotH265ProfileChanged(ui->cmbProfileH265->currentIndex());

    connect(ui->btnHdrMetadata, SIGNAL(clicked()), SLOT(slotEditHDRMetadata()));

    connect(ui->cmbPaletteuseDitherGIF, 
            SIGNAL(currentIndexChanged(int)),
            SLOT(slotBayerFilterSelected(int)));
    
    slotBayerFilterSelected(ui->cmbPaletteuseDitherGIF->currentIndex());
    
    setSupportsHDR(false);
}

KisVideoExportOptionsDialog::~KisVideoExportOptionsDialog()
{
    delete ui;
}

void KisVideoExportOptionsDialog::setSupportsHDR(bool value)
{
    m_d->supportsHDR = value;
    slotH265ProfileChanged(ui->cmbProfileH265->currentIndex());
}

KisPropertiesConfigurationSP KisVideoExportOptionsDialog::configuration() const
{
    KisPropertiesConfigurationSP cfg(new KisPropertiesConfiguration());

    cfg->setProperty("CodecId", currentCodecId());
    cfg->setProperty("h264PresetIndex", ui->cmbPresetH264->currentIndex());
    cfg->setProperty("h264ConstantRateFactor", ui->intCRFH264->value());
    cfg->setProperty("h264ProfileIndex", ui->cmbProfileH264->currentIndex());
    cfg->setProperty("h264TuneIndex", ui->cmbTuneH264->currentIndex());

    cfg->setProperty("h265PresetIndex", ui->cmbPresetH265->currentIndex());
    cfg->setProperty("h265ConstantRateFactor", ui->intCRFH265->value());
    cfg->setProperty("h265ProfileIndex", ui->cmbProfileH265->currentIndex());
    cfg->setProperty("h265TuneIndex", ui->cmbTuneH265->currentIndex());
    cfg->setProperty("h265UseHDRMetadata", ui->chkUseHDRMetadata->isChecked());

    cfg->setProperty("TheoraBitrate", ui->intBitrate->value());
    cfg->setProperty("CustomLineValue", ui->txtCustomLine->text());
    cfg->setProperty("customUserOptions", customUserOptions().join(' '));

    cfg->setPrefixedProperties("hdrMetadata/", m_d->hdrMetadataOptions.toProperties());

    return cfg;
}

KisVideoExportOptionsDialog::ContainerType KisVideoExportOptionsDialog::mimeToContainer(const QString &mimeType)
{
    if (mimeType == "video/webm") {
        return ContainerType::WEBM;
    } else if (mimeType == "video/ogg") {
        return ContainerType::OGV;
    } else if (mimeType == "image/gif") {
        return ContainerType::GIF;
    } else if (mimeType == "image/apng") {
        return ContainerType::APNG;
    } else if (mimeType == "image/webp") {
        return ContainerType::WEBP;
    }

    return ContainerType::DEFAULT;
}

void KisVideoExportOptionsDialog::slotCustomLineToggled(bool value)
{
    QString customLine = m_d->currentCustomLine;

    if (m_d->currentCustomLine.isEmpty() && value) {
        customLine = generateCustomLine().join(" ");
    } else if (!value) {
        customLine = QString();
        m_d->currentCustomLine = QString();
    }

    ui->txtCustomLine->setText(customLine);

    ui->stackedWidget->setEnabled(!value);
    ui->txtCustomLine->setEnabled(value);
    ui->btnResetCustomLine->setEnabled(value);
}

void KisVideoExportOptionsDialog::slotResetCustomLine()
{
    ui->txtCustomLine->setText(generateCustomLine().join(" "));
    slotSaveCustomLine();
}

void KisVideoExportOptionsDialog::slotCodecSelected(int index)
{
    const QString codec = m_d->codecs[index].id();

    if (codec == "libx264") {
        ui->stackedWidget->setCurrentIndex(CODEC_H264);
    } else if (codec == "libx265") {
        ui->stackedWidget->setCurrentIndex(CODEC_H265);
    } else if (codec == "libtheora") {
        ui->stackedWidget->setCurrentIndex(CODEC_THEORA);
    } else if (codec == "libvpx") {
        ui->stackedWidget->setCurrentIndex(CODEC_VP9);
    } else if (codec == "gif") {
        ui->stackedWidget->setCurrentIndex(CODEC_GIF);
    } else if (codec == "apng") {
        ui->stackedWidget->setCurrentIndex(CODEC_APNG);
    } else if (codec == "libwebp") {
        ui->stackedWidget->setCurrentIndex(CODEC_WEBP);
    }
}

void KisVideoExportOptionsDialog::slotSaveCustomLine()
{
    m_d->currentCustomLine = ui->txtCustomLine->text();
}

QStringList KisVideoExportOptionsDialog::customUserOptions() const
{
    return ui->chkCustomLine->isChecked() ?
        ui->txtCustomLine->text().split(" ", QString::SkipEmptyParts) :
                generateCustomLine();
}

QString KisVideoExportOptionsDialog::customUserOptionsString() const
{
    return customUserOptions().join(' ');
}

bool KisVideoExportOptionsDialog::videoConfiguredForHDR() const
{
    return currentCodecId() == "libx265" &&
        ui->chkUseHDRMetadata->isEnabled() &&
            ui->chkUseHDRMetadata->isChecked();
}

void KisVideoExportOptionsDialog::setHDRConfiguration(bool value) {
    if (value && currentCodecId() != "libx265") {
        ui->cmbCodec->setCurrentIndex(m_d->codecs.indexOf(KoID("libx265")));
        ui->chkUseHDRMetadata->setEnabled(true);
    }

    //If HDR is enabled && the codec id is correct, we need to use main10.
    if (value && currentCodecId() == "libx265") {
        ui->cmbProfileH265->setCurrentIndex(m_d->profilesH265.indexOf(KoID("main10")));
    }

    ui->chkUseHDRMetadata->setChecked(value);
}

int findIndexById(const QString &id, const QVector<KoID> &ids)
{
    int index = -1;
    auto it = std::find_if(ids.begin(), ids.end(),
                           [id] (const KoID &item) { return item.id() == id; });
    if (it != ids.end()) {
        index = std::distance(ids.begin(), it);
    }

    return index;
}

void KisVideoExportOptionsDialog::setConfiguration(const KisPropertiesConfigurationSP cfg)
{
    ui->cmbPresetH264->setCurrentIndex(cfg->getInt("h264PresetIndex", 5));
    ui->intCRFH264->setValue(cfg->getInt("h264ConstantRateFactor", 23));
    ui->cmbProfileH264->setCurrentIndex(cfg->getInt("h264ProfileIndex", 0));
    ui->cmbTuneH264->setCurrentIndex(cfg->getInt("h264TuneIndex", 1));

    ui->cmbPresetH265->setCurrentIndex(cfg->getInt("h265PresetIndex", 5));
    ui->intCRFH265->setValue(cfg->getInt("h265ConstantRateFactor", 23));
    ui->cmbProfileH265->setCurrentIndex(cfg->getInt("h265ProfileIndex", 0));
    ui->cmbTuneH265->setCurrentIndex(cfg->getInt("h265TuneIndex", 1));
    ui->chkUseHDRMetadata->setChecked(cfg->getBool("h265UseHDRMetadata", false));

    ui->intBitrate->setValue(cfg->getInt("TheoraBitrate", 5000));

    m_d->currentCustomLine = cfg->getString("CustomLineValue", QString());
    ui->chkCustomLine->setChecked(!m_d->currentCustomLine.isEmpty());
    slotCustomLineToggled(ui->chkCustomLine->isChecked());

    const QString codecId = cfg->getString("CodecId", "");

    const int index = qMax(0, findIndexById(codecId, m_d->codecs));
    ui->cmbCodec->setCurrentIndex(index);
    slotCodecSelected(index);

    slotH265ProfileChanged(ui->cmbProfileH265->currentIndex());

    KisPropertiesConfigurationSP metadataProperties = new KisPropertiesConfiguration();
    cfg->getPrefixedProperties("hdrMetadata/", metadataProperties);
    m_d->hdrMetadataOptions.fromProperties(metadataProperties);
}

QStringList KisVideoExportOptionsDialog::generateCustomLine() const
{
    QStringList options;

    if (currentCodecId() == "libx264") {
        options << "-crf" << QString::number(ui->intCRFH264->value());

        const int presetIndex = ui->cmbPresetH264->currentIndex();
        options << "-preset" << m_d->presets[presetIndex].id();

        const int profileIndex = ui->cmbProfileH264->currentIndex();
        options << "-profile:v" << m_d->profilesH264[profileIndex].id();

        if (m_d->profilesH264[profileIndex].id() == "high422") {
            options << "-pix_fmt" << "yuv422p";
        } else if (m_d->profilesH264[profileIndex].id() == "high444") {
            options << "-pix_fmt" << "yuv444p";
        } else {
            options << "-pix_fmt" << "yuv420p";
        }

        // Disabled! see the comment in c-tor!
        //const int tuneIndex = ui->cmbTune->currentIndex();
        //options << "-tune" << m_d->tunes[tuneIndex].id();

    } else if (currentCodecId() == "libx265") {
        const bool enableHDR =
            ui->chkUseHDRMetadata->isEnabled() &&
            ui->chkUseHDRMetadata->isChecked();

        if (enableHDR) {
            options << "-colorspace" << "bt2020c"
                    << "-color_trc" << "smpte2084"
                    << "-color_primaries" << "bt2020";
        }

        options << "-c:v" << "libx265";
        options << "-crf" << QString::number(ui->intCRFH265->value());

        const int presetIndex = ui->cmbPresetH265->currentIndex();
        options << "-preset" << m_d->presets[presetIndex].id();

        const int profileIndex = ui->cmbProfileH265->currentIndex();
        options << "-profile:v" << m_d->profilesH265[profileIndex].id();

        if (m_d->profilesH265[profileIndex].id() == "main") {
            options << "-pix_fmt" << "yuv420p";
        } else if (m_d->profilesH265[profileIndex].id() == "main10") {
            options << "-pix_fmt" << "yuv420p10le";
        } else {
            KIS_SAFE_ASSERT_RECOVER_NOOP(0 && "Unknown profile selected for h265 encoder");
        }

        if (enableHDR) {
            const QString metadataLine = m_d->hdrMetadataOptions.generateFFMpegOptions();
            options << metadataLine.split(" ");
        }

    } else if (currentCodecId() == "libtheora") {
        options << "-b" << QString::number(ui->intBitrate->value()) + "k";
    } else if (currentCodecId() == "libvpx") {
        options << "-vcodec" << "libvpx-vp9";
        if (ui->vp9Lossless->isChecked()) {
            options << "-lossless" <<  "1";
        } else {
            options << "-b:v" << QString::number(ui->vp9Mbits->value()) + "M";
        }
    } else if (currentCodecId() == "gif") {
        const QString ditherFilterString = m_d->paletteuseDither[ ui->cmbPaletteuseDitherGIF->currentIndex() ].id();

        options << "-f" << "gif"
                << "-loop" << ( ui->gifLoop->isChecked() ? "0":"-1" )
                << "-gifflags" << ( ui->gifTransDiff->isChecked() ? "+transdiff":"-transdiff" )
                << "-palettegen" << QString("palettegen=stats_mode=%1%2") 
                                            .arg(m_d->palettegenStatsMode[ ui->cmbPalettegenStatsModeGIF->currentIndex() ].id())
                                            .arg(":reserve_transparent=" + QString(ui->gifReserveTransparent->isChecked() ? "1":"0"))
                << "-lavfi" << QString("[0:v][1:v]paletteuse=dither=%1%2%3")
                                            .arg(ditherFilterString)
                                            .arg(ditherFilterString == "bayer" ? (QString(":bayer_scale=%1").arg(ui->intPaletteuseBayerScaleGIF->value()) ):"" )
                                            .arg(":diff_mode=" + m_d->paletteuseDiffMode[ ui->cmbPaletteuseDiffModeGIF->currentIndex() ].id() );
                
    } else if (currentCodecId() == "apng") {
        const int predIndex = ui->cmbPredAPNG->currentIndex();
        
        options << "-f" << "apng"
                << "-pred" << m_d->predAPNG[predIndex].id()
                << "-plays" << ( ui->apngLoop->isChecked() ? "0":"1" );
                
    } else if (currentCodecId() == "libwebp") {
        const int presetIndex = ui->cmbPresetWEBP->currentIndex();
        
        options << "-f" << "webp"
                << "-lossless" << ( ui->webpLossless->isChecked() ? "1":"0" )
                << "-compression_level" << QString::number(ui->intCompressWEBP->value())
                << "-q:v" << QString::number(ui->intQscaleWEBP->value())
                << "-preset" << m_d->presetWEBP[presetIndex].id()
                << "-loop" << ( ui->webpLoop->isChecked() ? "0":"1" );
        
    }

    return options;
}

QString KisVideoExportOptionsDialog::currentCodecId() const
{
    return m_d->codecs[ui->cmbCodec->currentIndex()].id();
}

void KisVideoExportOptionsDialog::slotH265ProfileChanged(int index)
{
    const bool enableHDR =
        m_d->supportsHDR &&
        index >= 0 &&
        m_d->profilesH265[index].id() == "main10";

    ui->chkUseHDRMetadata->setEnabled(enableHDR);
    ui->btnHdrMetadata->setEnabled(enableHDR && ui->chkUseHDRMetadata->isChecked());

    QString hdrToolTip;

    if (!m_d->supportsHDR) {
        hdrToolTip = i18nc("@info:tooltip", "Exported animation format does not support HDR");
    } else if (!enableHDR) {
        hdrToolTip = i18nc("@info:tooltip", "HDR metadata available only with \"main10\" profile");
    }

    ui->chkUseHDRMetadata->setToolTip(hdrToolTip);
    ui->btnHdrMetadata->setToolTip(hdrToolTip);
}

void KisVideoExportOptionsDialog::slotEditHDRMetadata()
{
    VideoHDRMetadataOptionsDialog dlg(this);
    dlg.setHDRMetadataOptions(m_d->hdrMetadataOptions);

    if (dlg.exec() == QDialog::Accepted) {
        m_d->hdrMetadataOptions = dlg.hdrMetadataOptions();
    }
}

void KisVideoExportOptionsDialog::slotBayerFilterSelected(int index)
{
    const bool enableBayer = m_d->paletteuseDither[ index ].id() == "bayer";
    ui->lblPaletteuseBayerScaleGIF->setEnabled( enableBayer );
    ui->intPaletteuseBayerScaleGIF->setEnabled( enableBayer );
}
