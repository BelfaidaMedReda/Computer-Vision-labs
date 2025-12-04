#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <random>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match {
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}


FMatrix<float,3,3> init_NormMatrix(float scale, float u_mean, float v_mean) {
    FMatrix<float,3,3> res = FMatrix<float,3,3>::Zero();
    res(0,0) = scale; res(1,1) = scale; res(2,2) = 1;
    res(0,2) = - u_mean * scale; res(1,2) = - v_mean * scale;
    return res;
}


// Implementing Hardley Normalization
// Ref. In the defense of Eight Points algorithm
FVector<FMatrix<float,3,3>,2> applyHartleyNormalization(vector<Match>& candidate_matches) {
    int n = candidate_matches.size();
    float u_mean = 0, v_mean = 0, v_mean_prime = 0, u_mean_prime = 0;
    float dist = 0, dist_prime = 0;
    auto U = Vector<float>::Zero(n);
    auto U_prime = Vector<float>::Zero(n);
    auto V = Vector<float>::Zero(n);
    auto V_prime = Vector<float>::Zero(n);
    for (int i =0; i < n; i++) {
        U[i] = candidate_matches[i].x1;
        U_prime[i] = candidate_matches[i].x2;
        V[i] = candidate_matches[i].y1;
        V_prime[i] = candidate_matches[i].y2;
    }
    u_mean = sum(U) / U.size(); v_mean = sum(V) / V.size();
    u_mean_prime = sum(U_prime) / U_prime.size(); v_mean_prime = sum(V_prime) / V_prime.size();
    for (int i =0; i < n; i++) {
        dist += sqrt(pow(U[i] - u_mean,2) + pow(V[i] - v_mean,2));
        dist_prime += sqrt(pow(U_prime[i] - u_mean_prime,2) + pow(V_prime[i] - v_mean_prime,2));
    }
    dist /= n;  dist_prime /= n;
    float s = sqrt(2) / dist;
    float s_prime = sqrt(2) / dist_prime;

    FMatrix<float,3,3> T = init_NormMatrix(s, u_mean, v_mean);
    FMatrix<float,3,3> T_prime = init_NormMatrix(s_prime,u_mean_prime,v_mean_prime);

    for (Match& m : candidate_matches) {
        FVector<float, 3> x = FVector<float, 3>(m.x1, m.y1, 1);
        FVector<float, 3> x_prime = FVector<float, 3>(m.x2, m.y2, 1);
        FVector<float, 3> x_norm = T * x;
        FVector<float, 3> x_prime_norm = T_prime * x_prime;
        m.x1 = x_norm.x(); m.y1 = x_norm.y(); m.x2 = x_prime_norm.x(); m.y2 = x_prime_norm.y();
    }
    return {T,T_prime};
}

FMatrix<float,3,3> estimateF(vector<Match> candidate_matches) {
    int n  = candidate_matches.size();
    float x1 = 0, x2 = 0, y1 = 0, y2 =0;
    int min_idx;
    FMatrix<float,3,3> candidateF;
    FVector<FMatrix<float,3,3>,2> NormMatrices;
    FMatrix<float,3,3> T, T_prime;
    // DO NOT FORGET NORMALIZATION OF POINTS
    // WE USE HARTLEY NORMALIZATION
    NormMatrices = applyHartleyNormalization(candidate_matches);
    T = NormMatrices[0];
    T_prime = NormMatrices[1];

    Matrix<float> A = Matrix<float>::Zero(n,9);
    for (int j=0; j < n; j++) {
        x1 = candidate_matches[j].x1;
        x2 = candidate_matches[j].x2;
        y2 = candidate_matches[j].y2;
        y1 = candidate_matches[j].y1;
        A(j,0) = x1 * x2;
        A(j,1) = x1 * y2;
        A(j,2) = x1;
        A(j,3) = x2 * y1;
        A(j,4) = y1 * y2;
        A(j,5) = y1;
        A(j,6) = x2;
        A(j,7) = y2;
        A(j,8) = 1;
    }

    // Apply SVD
    Matrix<float> U;
    Matrix<float> Vt;
    Vector<float> S;
    svd(A,U,S,Vt, true);
    min_idx = static_cast<int>(std::distance(std::begin(S), std::min_element(std::begin(S), std::end(S))));
    Vector<float> res = Vt.getRow(min_idx);

    // Fundamental Matrix Initialization
    candidateF(0,0) = res[0];
    candidateF(0,1) = res[1];
    candidateF(0,2) = res[2];
    candidateF(1,0) = res[3];
    candidateF(1,1) = res[4];
    candidateF(1,2) = res[5];
    candidateF(2,0) = res[6];
    candidateF(2,1) = res[7];
    candidateF(2,2) = res[8];

    // Apply SVD to F to force rank 2 constraint

    FMatrix<float,3,3> new_U;
    FMatrix<float,3,3> new_Vt;
    FMatrix<float,3,3> Sigma = FMatrix<float,3,3>::Zero();
    FVector<float,3> new_S;

    svd(candidateF,new_U,new_S,new_Vt, true);
    min_idx = static_cast<int>(std::distance(std::begin(new_S), std::min_element(std::begin(new_S), std::end(new_S))));
    new_S[min_idx] = 0;

    for (int k = 0; k < 3; k++) {
        Sigma(k,k) = new_S[k];
    }
    candidateF = new_U * Sigma * new_Vt;

    // Denormalization of F
    candidateF = transpose(T) * candidateF * T_prime;

    return candidateF;

}

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {
    vector<Match> all= matches;
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter=100000; // Adjusted dynamically
    int N_necessary = Niter; // Dynamic Niter
    float w_hat; // Useful to compute propotion of inliers
    FMatrix<float,3,3> bestF;
    vector<int> bestInliers;

    for (int i=0; i< Niter; i++) {
        vector<int> currentInliers;
        FMatrix<float,3,3> candidateF;
        vector<Match> candidate_matches;
        int random_idx;
        float x1 = 0, x2 = 0, y1 = 0, y2 = 0;
        float dist;
        Match newMatch;

        // Sampling process
        for (int j=0; j < 8; j++) {
            random_idx = rand() % matches.size();
            newMatch = matches[random_idx];
            candidate_matches.push_back(newMatch);
        }

        // Estimation process
        candidateF = estimateF(candidate_matches);

        // RANSAC Discrimination Logic
        for (int j = 0; j < matches.size(); j++) {
            x1 = matches[j].x1;
            x2 = matches[j].x2;
            y1 = matches[j].y1;
            y2 = matches[j].y2;
            FVector<float,3> x(x1,y1,1);
            FVector<float,3> x_prime(x2,y2,1);
            FVector<float,3> l_prime = transpose(candidateF) * x;
            dist = abs(sum(mult(x_prime,l_prime)) / norm2(l_prime));
            if (dist <= distMax) {
                currentInliers.push_back(j);
            }
        }

        if (currentInliers.size() > bestInliers.size()) {
            bestInliers.clear();
            bestInliers.swap(currentInliers);
            w_hat = bestInliers.size() / all.size();
            N_necessary = static_cast<int>(log(BETA) / log(1 - pow(w_hat,8)));
            matches.clear();
            for(size_t k=0; k < bestInliers.size(); k++)
                matches.push_back(all[bestInliers[k]]);
            bestF = estimateF(matches);
        }

        if (i >= N_necessary) {
            break;
        }
    }
    return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {
    const int w1 = I1.width();
    const int w2 = I2.width();
    while(true) {
        int x,y;
        if (getMouse(x,y) == 1) {
            fillCircle(x+0, y, 2, RED);
            bool OnLeft = x <= w1 - 1;
            auto new_x = static_cast<float>(x);
            auto new_y = static_cast<float>(y);
            if (OnLeft) {
                FVector<float,3> X(new_x, new_y,1);
                FVector<float,3> l_prime = transpose(F) * X;
                float y1 = - l_prime[2] / l_prime[1];
                float y2 = - (l_prime[2] + l_prime[0] * static_cast<float>(w2)) / l_prime[1];
                drawLine(w1, static_cast<int>(y1), w2 + w1,static_cast<int>(y2), RED);
            }
            else {
                FVector<float,3> X_prime(new_x - static_cast<float>(w1), new_y,1);
                FVector<float,3> l = F * X_prime;
                float y1 = - l[2] / l[1];
                float y2 = - (l[2] + l[0] * static_cast<float>(w1)) / l[1];
                drawLine(0, static_cast<int>(y1), w1,static_cast<int>(y2), RED);
            }
        }
        if(getMouse(x,y) == 3)
            break;

    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    string s1 = argc>1? argv[1]: srcPath("im1.jpg");
    string s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    const int n = (int)matches.size();
    cout << " matches: " << n << endl;
    drawString(100,20,std::to_string(n)+ " matches",RED);
    click();
    
    FMatrix<float,3,3> F = computeF(matches);
    cout << "F="<< endl << F;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c(rand()%256,rand()%256,rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);        
    }
    drawString(100, 20, to_string(matches.size())+"/"+to_string(n)+" inliers", RED);
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}
