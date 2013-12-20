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
    double h = 0.001;
    double damping = 0.01;
    double crushspeed = 20;
    double crushFraction = 0.10;
    double stretchStiffness = 10000000;
    double bendStiffness    = 10000;

    VectorXd q(numdofs());
    VectorXd g(numedges());
    dofsFromGeometry(q, g);
    VectorXd undefq = q;


    VectorXd v(numdofs());
    v.setZero();
    ShellForces sf(*mesh_, undefq, stretchStiffness, bendStiffness);
    SparseMatrix<double> Minv;
    SparseMatrix<double> M;
    buildMassMatrix(g, M);
    buildInvMassMatrix(g, Minv);


    // find top and bottom
    double miny = std::numeric_limits<double>::infinity();
    double maxy = -std::numeric_limits<double>::infinity();
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        double y = mesh_->point(vi.handle())[1];
        miny = std::min(miny, y);
        maxy = std::max(maxy, y);
    }

    for(int i=0; h*i*crushspeed < crushFraction*(maxy-miny); i++)
    {
        q += h*v;

        // fix bottom and top
        for(int j=0; j<(int)mesh_->n_vertices(); j++)
        {
            double height = q[3*j+1];
            if(height < miny)
            {
                q[3*j+1] = miny;
                //DOFs on the floor become pinned there
                for(int k=0; k<3; k++)
                {
                    v[3*j+k] = 0;
                    Minv.coeffRef(3*j+k,3*j+k) = 0;
                }
            }
            if(height > maxy-h*i*crushspeed)
            {
                q[3*j+1] = maxy - h*i*crushspeed;
                v[3*j+1] = -crushspeed;
            }
        }


        v *= (1.0 - damping);

        VectorXd F;
        sf.getForce(q, F);
        v += h*Minv*F;
        dofsToGeometry(q, g);
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

