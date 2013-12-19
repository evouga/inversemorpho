#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include <string>
#include "mesh.h"

class Controller;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void setController(Controller &cont);

public slots:
    void setParameters(ProblemParameters params);
    void showError(std::string error);
    void centerCamera(Eigen::Vector3d centroid, double radius);
    void repaintMesh();
    void tick();

private slots:
    void on_actionExit_triggered();

    void on_actionReset_Camera_triggered();

    void on_actionTake_Screenshot_triggered();

    void on_wireframeCheckBox_clicked();

    void on_smoothShadeCheckBox_clicked();

    void on_actionExport_OBJ_triggered();

    void on_actionImport_OBJ_triggered();
    void on_relaxEmbeddingButton_clicked();

private:
    void updateGL();
    void saveScreenshot();
    void saveScreenshot(const std::string &filename);
    ProblemParameters getParameters();
    std::string launchImportOBJDialog();

    Ui::MainWindow *ui;
    Controller *cont_;
    QTimer *repainttimer_;
};

#endif // MAINWINDOW_H
