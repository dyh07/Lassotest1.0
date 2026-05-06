#include "mainwindow.h"
#include <QWebEngineSettings>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    // 创建网页控件
    webView = new QWebEngineView(this);
    setCentralWidget(webView);
    resize(2000, 1000);

    // 仅开启基础JS
    webView->settings()->setAttribute(QWebEngineSettings::JavascriptEnabled, true);

    // 加载Python网页
    webView->setUrl(QUrl("http://127.0.0.1:15114"));
}

MainWindow::~MainWindow()
{
    delete webView;
}