#include "mesh.h"
#include <iomanip>
#include "controller.h"
#include <Eigen/Dense>
#include "shellforces.h"
#include <QThread>

using namespace std;
using namespace Eigen;
using namespace OpenMesh;

bool Mesh::simulate(Controller &cont)
{
    VectorXd q(numdofs());
    VectorXd g(numedges());
    dofsFromGeometry(q, g);

    double h = 0.01;
    VectorXd v(numdofs());
    v.setZero();
    int numsteps = 10000;
    ShellForces sf(*mesh_, q, 100, 1);
    SparseMatrix<double> Minv;
    SparseMatrix<double> M;
    buildMassMatrix(g, M);
    buildInvMassMatrix(g, Minv);

    double damping = 0.1;

    // find top vertices
    double miny = std::numeric_limits<double>::infinity();
    double maxy = -std::numeric_limits<double>::infinity();
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        double y = mesh_->point(vi.handle())[1];
        miny = std::min(miny, y);
        maxy = std::max(maxy, y);
    }

    vector<int> topverts;
    vector<int> bottomverts;
    double tol = 1e-4;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        double y = mesh_->point(vi.handle())[1];
        if(fabs(y-miny) < tol)
            bottomverts.push_back(vi.handle().idx());
        if(fabs(y-maxy) < tol)
            topverts.push_back(vi.handle().idx());
    }

    Vector3d gravity(0,-20,0);
    VectorXd extForce(q.size());
    extForce.setZero();
    for(int i=0; i<(int)topverts.size(); i++)
    {
        extForce.segment<3>(3*topverts[i]) = gravity;
    }
    extForce = M*extForce;

    for(int i=0; i<(int)bottomverts.size(); i++)
    {
        for(int j=0; j<3; j++)
            Minv.coeffRef(3*bottomverts[i]+j, 3*bottomverts[i]+j) = 0;
    }

    for(int i=0; i<numsteps; i++)
    {
        q += h*v;

        v *= (1.0 - damping);

        VectorXd F;
        sf.getForce(q, F);
        F += extForce;
        v += h*Minv*F;
        dofsToGeometry(q, g);
        QThread::msleep(100);
    }

    cont.updateGL();
    return true;
}

void Mesh::buildInvMassMatrix(const VectorXd &g, Eigen::SparseMatrix<double> &Minv) const
{
    Minv.resize(numdofs(), numdofs());
    vector<Tr> entries;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        int vidx = vi.handle().idx();
        double area = barycentricDualArea(g, vidx);
        double invmass = 1.0/area;///params_.rho/params_.h/params_.scale/params_.scale/params_.scale;
//        if(mesh_->is_boundary(vi.handle()))
//            invmass = 0;
        for(int i=0; i<3; i++)
            entries.push_back(Tr(3*vidx+i, 3*vidx+i, invmass));
    }

    Minv.setFromTriplets(entries.begin(), entries.end());
}

void Mesh::buildMassMatrix(const VectorXd &g, Eigen::SparseMatrix<double> &Minv) const
{
    Minv.resize(numdofs(), numdofs());
    vector<Tr> entries;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        int vidx = vi.handle().idx();
        double area = barycentricDualArea(g, vidx);
        double invmass = area;///params_.rho/params_.h/params_.scale/params_.scale/params_.scale;
        for(int i=0; i<3; i++)
            entries.push_back(Tr(3*vidx+i, 3*vidx+i, invmass));
    }

    Minv.setFromTriplets(entries.begin(), entries.end());
}

double Mesh::barycentricDualArea(const VectorXd &g, int vidx) const
{
    double result = 0;
    OMMesh::VertexHandle vh = mesh_->vertex_handle(vidx);
    for(OMMesh::VertexFaceIter vfi = mesh_->vf_iter(vh); vfi; ++vfi)
    {
        result += restFaceArea(g, vfi.handle().idx());
    }
    return result/3.0;
}

double Mesh::restFaceArea(const VectorXd &g, int fidx) const
{
    OMMesh::FaceHandle fh = mesh_->face_handle(fidx);
    double gs[3];
    int idx=0;
    for(OMMesh::FaceEdgeIter fei = mesh_->fe_iter(fh); fei; ++fei)
    {
        gs[idx++] = g[fei.handle().idx()];
    }

    double s = 0.5*(gs[0]+gs[1]+gs[2]);
    return sqrt(s*(s-gs[0])*(s-gs[1])*(s-gs[2]));
}



double Mesh::faceArea(const VectorXd &q, int fidx) const
{
    FaceHandle fh = mesh_->face_handle(fidx);
    int verts[3];
    int idx=0;
    for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(fh); fvi; ++fvi)
    {
        verts[idx++] = fvi.handle().idx();
    }

    Vector3d q0 = q.segment<3>(3*verts[0]);
    Vector3d q1 = q.segment<3>(3*verts[1]);
    Vector3d q2 = q.segment<3>(3*verts[2]);

    double A = ((q1-q0).cross(q2-q0)).norm();
    return 0.5*A;
}

