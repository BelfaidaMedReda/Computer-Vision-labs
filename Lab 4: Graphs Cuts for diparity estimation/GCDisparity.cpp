// Imagine++ project
// Project:  GraphCutsDisparity
// Author:   Renaud Marlet/Pascal Monasse
// Student: TODO

#include <Imagine/Images.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>
#include "maxflow/graph.h"
using namespace Imagine;
using namespace std;

typedef Image<octet> octetImage;
typedef Image<double> doubleImage;

// Default data
const string DEF_im1=srcPath("assets/face0.png"), DEF_im2=srcPath("assets/face1.png");
static int dmin=10, dmax=55; // Min and max disparities

// Parameters of the algorithm
// OPTIMIZATION: to make the program faster, a zoom factor is used to
// down-sample the input images on the fly. You will
// only look at pixels (win+zoom*i,win+zoom*j) with win the radius of patch.
const int win = (7-1)/2;    // Correlation patches of size (2n+1)*(2n+1)
const float lambdaf = 0.1;  // Weight of regularization (smoothing) term
const int zoom = 2;         // Zoom factor (to speedup computations)
const float sigma = 6/zoom; // Gaussian blur parameter for disparity
// Energy discretization precision (as we build a graph with 'int' weights)
const int wcc = std::max(1+int(1/lambdaf),20);
const int lambda = lambdaf * wcc; // Regularization term (must be >= 1)

Image<Color> displayDisp(const Image<int>& disp) {
    Image<Color> im(disp.width(), disp.height());
    for(int j=0; j<disp.height(); j++)
        for(int i=0; i<disp.width(); i++) {
            if(disp(i,j)<dmin || disp(i,j)>dmax)
                im(i,j) = CYAN;
            else {
                int g = 255*(disp(i,j)-dmin)/(dmax-dmin);
                im(i,j)= Color(g,g,g);
            }
        }
    return im;
}

#ifdef IMAGINE_OPENGL
void doc() {
    cout << "***** 3D mesh renderings *****" << endl;
    cout << "- Button 1: toggle textured or gray rendering" << endl;
    cout << "- SHIFT+Button 1: rotate" << endl;
    cout << "- SHIFT+Button 3: translate" << endl;
    cout << "- Mouse wheel: zoom" << endl;
    cout << "- SHIFT+a: zoom out" << endl;
    cout << "- SHIFT+z: zoom in" << endl;
    cout << "- SHIFT+r: recenter camera" << endl;
    cout << "- SHIFT+m: toggle solid/wire/points mode" << endl;
    cout << "- Button 3: exit" << endl;
}

// Display 3D mesh renderings
void control3D(const Mesh& Mt, const Mesh& Mg) {
    setActiveWindow( openWindow3D(512, 512, "3D") );
    showMesh(Mt);
    bool textured = true;
    while (true) {
        Event evt;
        getEvent(5, evt);
        // On mouse button 1
        if (evt.type == EVT_BUT_ON && evt.button == 1) {
            // Toggle textured rendering and gray rendering
            if (textured) {
                hideMesh(Mt,false);
                showMesh(Mg,false);
            } else {
                hideMesh(Mg,false);
                showMesh(Mt,false);
            }
            textured = !textured;
        }
        // On mouse button 3
        if (evt.type == EVT_BUT_ON && evt.button == 3)
            break;
    }
}
#endif

void show3D(const octetImage& I, const doubleImage& D, int zoom) {
#ifdef IMAGINE_OPENGL
    cout << "Click to compute depth map and 3D mesh renderings... " << flush;
    click();

    // Compute 3D point cloud: magic constants depending on camera pose
    const float f  = 750;  // Focal
    const float d0 = 100;  // Disparity for infinite depth (due to crop)
    const float B  = -0.20;// Baseline

    const int nx=D.width(), ny=D.height();
    FMatrix<float,3,3> K(0.0f);
    K(0,0)= -f/zoom; K(0,2)=nx/2;
    K(1,1)=  f/zoom; K(1,2)=ny/2;
    K(2,2)=1.0f;
    K = inverse(K);
    K /= K(2,2);
    Array<FloatPoint3> p(nx*ny);
    Array<Color> pcol(nx*ny);
    for(int j=0; j<ny; j++)
        for(int i=0; i<nx; i++) {
            float z = f*B/(d0+D(i,j));
            FloatPoint3 pt((float)i,(float)j,1.0f);
            p[i+nx*j] = K*pt*z;
            pcol[i+nx*j] = Color(I(i*zoom,j*zoom));
        }
    // Create mesh from 3D point cloud
    Array<Triangle> t(2*(nx-1)*(ny-1));
    Array<Color> tcol(2*(nx-1)*(ny-1));
    for(int j=0; j+1<ny; j++)
        for(int i=0; i+1<nx; i++) {
            // Create triangles with next pixels in line/column
            t[2*(i+j*(nx-1))]   = Triangle(i+nx*j,   i+1+nx*j,     i+nx*(j+1));
            t[2*(i+j*(nx-1))+1] = Triangle(i+1+nx*j, i+1+nx*(j+1), i+nx*(j+1));
            tcol[2*(i+j*(nx-1))]   = pcol[i+nx*j];
            tcol[2*(i+j*(nx-1))+1] = pcol[i+nx*j];
        }
    // Mesh with texture from original image
    Mesh Mt(p.data(), nx*ny, t.data(), 2*(nx-1)*(ny-1), 0, 0, FACE_COLOR);
    Mt.setColors(TRIANGLE, tcol.data());
    // Mesh with artificial light
    Mesh Mg(p.data(), nx*ny, t.data(), 2*(nx-1)*(ny-1), 0, 0,
            CONSTANT_COLOR, SMOOTH_SHADING);
    cout << "done" << endl;

    doc();
    control3D(Mt, Mg);
#else
    cout << "No 3D: Imagine++ not built with OpenGL support" << endl;
#endif
}

// Return image of mean intensity value over (2win+1)x(2win+1) patch
doubleImage meanImage(const octetImage& I) {
    int w = I.width(), h = I.height();
    doubleImage IM(w,h);
    double area = (2*win+1)*(2*win+1);
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            if (j-win<0 || j+win>=h || i-win<0 ||i+win>=w) {
                IM(i,j)=0;
                continue;
            }
            double sum = 0;
            for(int y=-win; y<=win; y++)
                for(int x=-win; x<=win; x++)
                    sum += I(x+i,y+j);
            IM(i,j) = sum / area;
        }
    return IM;
}

// Compute correlation between two pixels in images 1 and 2
double correl(const octetImage& I1,    // Image 1
              const doubleImage& I1M, // Image of mean value over patch
              const octetImage& I2,    // Image2
              const doubleImage& I2M, // Image of mean value over patch
              int u1, int v1,         // Pixel of interest in image 1
              int u2, int v2) {       // Pixel of interest in image 2
    double c=0;
    for(int y=-win; y<=win; y++)
        for(int x=-win; x<=win; x++) 
            c += (I1(u1+x,v1+y)-I1M(u1,v1)) * (I2(u2+x,v2+y)-I2M(u2,v2));
    return c / ((2*win+1)*(2*win+1));
}

// Compute ZNCC between two patches in images 1 and 2
double zncc(const octetImage& I1,    // Image 1
            const doubleImage& I1M, // Image of mean intensity value over patch
            const octetImage& I2,    // Image2
            const doubleImage& I2M, // Image of mean intensity value over patch
            int u1, int v1,         // Pixel of interest in image 1
            int u2, int v2) {       // Pixel of interest in image 2
    double var1 = correl(I1, I1M, I1, I1M, u1, v1, u1, v1);
    if(var1 == 0)
        return 0;
    double var2 = correl(I2, I2M, I2, I2M, u2, v2, u2, v2);
    if(var2 == 0)
        return 0;
    return correl(I1, I1M, I2, I2M, u1, v1, u2, v2) / sqrt(var1 * var2);
}

/// Encodes a coordinate (x, y, d) into a single index
int encode_coordinates(int x, int y, int d, int nx, int ny) {
    return d * nx * ny + y * nx + x;
}

double rho(double x) {
    if (x < 0) {
        return 1;
    }
    return sqrt(1 - x);
}

bool valid_patch(int u, int v, int w, int h) {
    return (u - win >= 0 &&
            u + win < w &&
            v - win >= 0 &&
            v + win < h);
}

/// Create graph
/// The graph library works with node numbers. To clarify the setting, create
/// a formula to associate a unique node number to a triplet (x,y,d) of pixel
/// coordinates and disparity.
/// The library assumes an edge consists of a pair of oriented edges, one in
/// each direction. Put correct weights to the edges, such as 0, INF, or
/// an intermediate weight.
void build_graph(Graph<int,int,int>& G,
                 const octetImage& I1, const octetImage& I2,
                 int nx, int ny, int nd) {
    const int INF= 100000000; // "Infinite" value for edge impossible to cut
    // Precompute images of mean intensity value over patch
    doubleImage I1M = meanImage(I1), I2M = meanImage(I2);
    int K = 1 + 4 * nd * lambda;
    (void)G.add_node(nd * nx * ny);
    // Creating horizental edges
    for(int d=0; d< nd - 1; d++) {
        int shifted_d = d + dmin + 1;
        for (int x=0; x<nx; x++) {
            for (int y=0; y<ny; y++) {
                int u = win + x * zoom;
                int v = win + y * zoom;
                int curr_pixel = encode_coordinates(x, y, d, nx, ny);
                int next_pixel = encode_coordinates(x, y, d + 1, nx, ny);
                int forward_weight;
                if(valid_patch(u,v,I1.width(),I1.height()) && valid_patch(u+shifted_d,v,I2.width(),I2.height())){
                    forward_weight =  K + wcc * rho(zncc(I1, I1M, I2, I2M, u, v, u + shifted_d, v));
                }
                else{
                    forward_weight = INF;
                }
                int backward_weight = 0;
                G.add_edge(curr_pixel, next_pixel, forward_weight, backward_weight);
            }
        }
    }
    // Creating n-links
    for (int d=0; d<nd; d++) {
        for (int x=0; x<nx; x++) {
            for (int y=0; y<(ny - 1); y++) {
                int curr_pixel = encode_coordinates(x, y, d, nx, ny);
                int next_pixel = encode_coordinates(x, y + 1, d, nx, ny);
                int forward_weight = lambda;
                int backward_weight = lambda;
                G.add_edge(curr_pixel, next_pixel, forward_weight, backward_weight);
            }
        }
        for (int x=0; x<(nx -1); x++) {
            for (int y=0; y<ny; y++) {
                int curr_pixel = encode_coordinates(x, y, d, nx, ny);
                int next_pixel = encode_coordinates(x + 1, y, d, nx, ny);
                int forward_weight = lambda;
                int backward_weight = lambda;
                G.add_edge(curr_pixel, next_pixel, forward_weight, backward_weight);
            }
        }
    }
    // Creating t-links
    for (int x=0; x<nx; x++) {
        for (int y=0; y<ny; y++) {
            int u = win + x * zoom;
            int v = win + y * zoom;
            int curr_pixel = encode_coordinates(x, y, 0, nx, ny);
            int source_weight;
            if(valid_patch(u,v,I1.width(),I1.height()) && valid_patch(u+ dmin,v,I2.width(),I2.height())){
                source_weight =   K + wcc * rho(zncc(I1, I1M, I2, I2M, u, v, u + dmin, v));
            }
            else{
                source_weight = INF;
            }
            int sink_weight = 0;
            G.add_tweights(curr_pixel, source_weight, sink_weight);
        }
    }
    for (int x=0; x<nx; x++) {
        for (int y=0; y<ny; y++) {
            int u = win + x * zoom;
            int v = win + y * zoom;
            int curr_pixel = encode_coordinates(x, y, nd - 1, nx, ny);
            int source_weight = 0;
            int sink_weight;
            if(valid_patch(u,v,I1.width(),I1.height()) && valid_patch(u + dmax,v,I2.width(),I2.height())){
                sink_weight =  K + wcc * rho(zncc(I1, I1M, I2, I2M, u, v, u + dmax, v));
            }
            else{
                sink_weight = INF;
            }
            G.add_tweights(curr_pixel, source_weight, sink_weight);
        }
    }
}

/// Extract disparity from minimum cut.
/// w2 is the width of image2, used to detect pixels with no valid patch.
doubleImage decode_graph(Graph<int,int,int>& G, int nx, int ny, int nd, int w2){
    int curr_pixel, next_pixel;
    int u, shifted_d;
    doubleImage D(nx,ny);
    int unassigned = 0;
    for (int x=0; x<nx; x++) {
        for (int y=0; y<ny; y++) {
            u = win + x * zoom;
            int first_pixel = encode_coordinates(x, y, 0, nx, ny);
            int last_pixel = encode_coordinates(x, y, nd - 1, nx, ny);
            if (G.what_segment(first_pixel) == Graph<int,int,int>::SINK){
                shifted_d = dmin + 0;
                D(x, y) = (u + shifted_d < w2)? dmin : 0;
            }
            else if (G.what_segment(last_pixel) == Graph<int,int,int>::SOURCE) {
                shifted_d = dmin + nd - 1;
                D(x, y) =  (u + shifted_d < w2)? shifted_d : 0;
            }
            else {
                for (int d=0; d< nd-1; d++) {
                    curr_pixel = encode_coordinates(x, y, d, nx, ny);
                    next_pixel =  encode_coordinates(x, y, d+1, nx, ny);
                    if (G.what_segment(curr_pixel) == Graph<int,int,int>::SOURCE && G.what_segment(next_pixel) == Graph<int,int,int>::SINK) {
                        shifted_d = dmin + d + 1;
                        D(x, y) =  (u + shifted_d < w2)? shifted_d : 0;
                        break;
                    }
                }
            }
            if (D(x,y) == 0) {
                unassigned++;
            }
        }
    }
    std::cout << "Unassigned pixels" << unassigned << std::endl;
    return D;
}
    
// Load two rectified images.
// Compute the disparity of image 2 w.r.t. image 1.
// Display disparity map.
// Display 3D mesh of corresponding depth map.
int main(int argc, char* argv[]) {
    if(argc!=1 && argc!=5) {
        cerr << "Usage: " << argv[0] << " im1 im2 dmin dmax" << endl;
        return 1;
    }
    string im1=DEF_im1, im2=DEF_im2;
    if(argc>1) {
        im1 = argv[1]; im2=argv[2]; dmin=stoi(argv[3]); dmax=stoi(argv[4]);
    }
    cout << "Loading images... " << flush;
    octetImage I1,I2;
    if(!load(I1, im1) || !load(I2,im2)) {
        cerr << "Error loading image files" << endl;
        return 1;
    }
    cout << "done" << endl;

    cout << "Parameters: " << "d=" << dmin << "..." << dmax
         << ", win="<<win << ", lambda="<<lambdaf << ", sigma="<<sigma
         << ", zoom="<<zoom << endl;

    cout << "Displaying images... " << flush;
    int w1=I1.width(), w2=I2.width(), h=I1.height();
    const int nd=dmax-dmin; // Disparity range
    openWindow(w1+w2, h);
    display(I1); display(I2,w1,0);
    cout << "done" << endl;

    // Zoomed image dim, disregarding borders (strips of width the patch radius)
    const int nx=(w1-2*win)/zoom, ny=(h-2*win)/zoom;

    cout << "Constructing graph (be patient)... " << flush;
    Graph<int,int,int> G(nx*ny*nd,2*nx*ny*nd);
    build_graph(G, I1, I2, nx, ny, nd);
    cout << "done" << endl;

    cout << "Computing minimum cut... " << flush;
    int f = G.maxflow();
    cout << "done" << endl << "  max flow = " << f << endl;

    cout << "Extracting disparity map from minimum cut... " << flush;
    doubleImage D=decode_graph(G, nx, ny, nd, I2.width());
    cout << "done" << endl;

    cout << "Displaying disparity map... " << flush;
    fillRect(0,0,w1,h,CYAN);
    display(enlarge(displayDisp(D),zoom),win,win);
    cout << "done" << endl;
    cout << "Click to compute and display blurred disparity map... " << flush;
    click();
    D=blur(D,sigma);
    display(enlarge(displayDisp(D),zoom),win,win);
    cout << "done" << endl;

    show3D(I1.getSubImage(win,win,w1-2*win,h-2*win), D, zoom);
    endGraphics();
    return 0;
}