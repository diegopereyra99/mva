// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>    // std::shuffle
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

FMatrix<float,3,3> FundamentalMatrix8(vector<Match> matches){

    FMatrix<float,9,9> U, Vt, A;
    A.fill(0);
    FVector<float,9> Sigmas;

    for(int i=0; i < 8; i++) {
        Match m = matches[i]; 

        A(i, 0) = m.x2 * m.x1;
        A(i, 1) = m.x2 * m.y1;
        A(i, 2) = m.x2;
        A(i, 3) = m.y2 * m.x1;
        A(i, 4) = m.y2 * m.y1;
        A(i, 5) = m.y2;
        A(i, 6) = m.x1;
        A(i, 7) = m.y1;
        A(i, 8) = 1;       
    }

    svd(A, U, Sigmas, Vt);

    // Last column of Vt is F    
    FVector<float,3> Sigmas_F(3);
    FMatrix<float,3,3> U_F, Vt_F, F;
    F(0,0) = Vt(8,0); 
    F(0,1) = Vt(8,1); 
    F(0,2) = Vt(8,2); 
    F(1,0) = Vt(8,3); 
    F(1,1) = Vt(8,4); 
    F(1,2) = Vt(8,5); 
    F(2,0) = Vt(8,6); 
    F(2,1) = Vt(8,7); 
    F(2,2) = Vt(8,8); 
    
    // SVD on F to then force the min s.v. to be 0
    svd(F, U_F, Sigmas_F, Vt_F);
    Sigmas_F[2] = 0;
    F = U_F * Diagonal(Sigmas_F) * Vt_F;

    return F;

}

vector<int> getInliers(FMatrix<float, 3, 3>& F, vector<Match>& matches, float distMax){
    FVector<float, 3> x, x_;
    vector<int> inliers;

    for (int i=0; i < matches.size(); i++) {
        x[0] = matches[i].x1;
        x[1] = matches[i].y1;
        x[2] = 1;
        x_[0] = matches[i].x2;
        x_[1] = matches[i].y2;
        x_[2] = 1;

        FVector<float, 3> Fx = F * x;
        float dist = abs(x_ * Fx);
        dist /= sqrt(Fx[0] * Fx[0] + Fx[1] * Fx[1]);
        if (dist <= distMax){
            inliers.push_back(i);
        }
    }

    return inliers;
}

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter=100000; // Adjusted dynamically
    FMatrix<float,3,3> bestF;
    vector<int> bestInliers;
    // --------------- TODO ------------
    // DO NOT FORGET NORMALIZATION OF POINTS

    const float NORM = 1e-3;
    const int N_matches =  matches.size();
    FMatrix<float, 3, 3> N;
    N.fill(0);
    N(0, 0) = NORM;
    N(1, 1) = NORM; 
    N(2, 2) = 1; 

    vector<Match> norm_matches;
    for(int i=0; i<N_matches; i++) {
        Match m = matches[i];
        m.x1 = m.x1 * NORM;
        m.y1 = m.y1 * NORM;
        m.x2 = m.x2 * NORM;
        m.y2 = m.y2 * NORM;
        norm_matches.push_back(m);
    }       

    // RANSAC iterations
    long int iter=0;
    while (iter < Niter){

        // Select 8 random indices for the matches
        vector<int> indices(8);
        vector<Match> random_matches(8);
        int c = 0;
        while (c < 8) {
            int ind = rand() % N_matches;
            if (find(indices.begin(), indices.end(), ind) == indices.end()){
                indices[c] = ind;
                random_matches[c] = norm_matches[ind];
                c++;
            }
        }
        
        // Calculate F
        FMatrix<float,3,3> F;
        F = FundamentalMatrix8(random_matches);
        F = N * F * N;

        // See how many inliers there are
        vector<int> iter_inliers;
        iter_inliers = getInliers(F, matches, distMax);

        // check if with the new F there are more inliers
        if (iter_inliers.size() > bestInliers.size()){
            cout << "iter=" << iter << " / num_inliers=" << iter_inliers.size() << endl;
            bestInliers = iter_inliers;
            bestF = F;
            float m = bestInliers.size();
            Niter = (int)(log(BETA) / log(1 - pow(m / N_matches, 8)));
        }

        iter++;
    }

    // Updating matches with inliers only
    vector<Match> all=matches;
    matches.clear();
    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);
    return bestF;
}





// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {
    while(true) {
        int x,y;
        if(getMouse(x,y) == 3)
            break;
        // --------------- TODO ------------

        IntPoint2 point{x, y};
        FVector<int, 3> X{x, y, 1};
        FVector<float, 3> epipolarLine;
        float a, b, x1, x2, y1, y2;
        Color color;
        
        if (x < I1.width()) {
            // Line in im2
            epipolarLine = F * X;
            x1 = 0;
            x2 = I2.width();
            color = RED;
            
        }
        else {
            // Line in im1
            X[0] -= I1.width();
            epipolarLine = transpose(F) * X;
            x1 = 0;
            x2 = I1.width();
            color = GREEN;
        }
    
        a = -epipolarLine[0] / epipolarLine[1];
        b = -epipolarLine[2] / epipolarLine[1];

        if (x < I1.width()){
            x1 += I1.width();
            x2 += I1.width();
        }
        y1 = b;
        y2 = a * x2 + b;

        drawLine(x1, y1, x2, y2, color);
        drawCircle(point, 3, color);
    }
}


int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

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
