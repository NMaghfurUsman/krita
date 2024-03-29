/*
 *  SPDX-FileCopyrightText: 2016 Boudewijn Rempt <boud@valdyas.org>
 *  SPDX-FileCopyrightText: 2021 L. E. Segovia <amy@amyspark.me>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 *
 */
#ifndef KISREMOTEFILEFETCHER_H
#define KISREMOTEFILEFETCHER_H

#include <QIODevice>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QObject>
#include <QUrl>

/**
 * @brief The KisRemoteFileFetcher class can fetch a remote file and blocks until the file is downloaded
 */
class KisRemoteFileFetcher : public QObject
{
    Q_OBJECT
public:
    explicit KisRemoteFileFetcher(QObject *parent = 0);
    ~KisRemoteFileFetcher() override;

    /// fetch the image. Shows a progress dialog
    bool fetchFile(const QUrl &remote, QIODevice *io);

    static QByteArray fetchFile(const QUrl &remote);


private Q_SLOTS:
    void error(QNetworkReply::NetworkError error);

private:
    QNetworkRequest *m_request;
    QNetworkReply *m_reply;
};

#endif // KISREMOTEFILEFETCHER_H
