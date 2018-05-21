#include "HelperFunctions.h"

namespace OPS{

SimulationState SimulationState::readFromFile(std::string file){
    size_t n, ns, s;
    double_t g, b;
    ifstream f( file );
    std::string line;
    // Read number of particles
    std::getline( f, line );
    std::getline( f, line );
    std::getline( f, line );
    n = std::stod(line);
    // Read nameSuffix
    std::getline( f, line );
    std::getline( f, line );
    ns = std::stod(line);
    // Read step number
    std::getline( f, line );
    std::getline( f, line );
    s = std::stod(line);
    // Read gamma or the FvK number
    std::getline( f, line );
    std::getline( f, line );
    g = std::stod(line);
    // Read beta or 1/temperature
    std::getline( f, line );
    std::getline( f, line );
    b = std::stod(line);
    // Read random engine state
    std::getline( f, line );
    std::getline( f, line );
    std::istringstream ss(line);
    std::mt19937 en;
    ss >> en;
    // Read random generator state
    std::getline( f, line );
    //std::getline( f, line );
    //std::istringstream ss(line);
    NormD rg;
    f >> rg;
    // Read neighbors
    vtkIdType id;
    IdList l(n);
    std::getline( f, line );
    std::getline( f, line );
    size_t counter = 0;
    while( counter < n &&  std::getline( f, line ) ){
        std::istringstream ss(line);
        for(auto i=0; i < numsPerLine; ++i){
            ss >> id;
            l[counter++] = id;
        }
    }
    // Read initial position of particles
    M3X ip(3,n);
    double_t xi;
    std::getline( f, line );
    counter = 0;
    while( counter < 3*n &&  std::getline( f, line ) ){
        std::istringstream ss(line);
        for(auto i=0; i < numsPerLine; ++i){
            ss >> xi;
            ip(counter++) = xi;
        }
    }
    // Read prevX vector
    Vec pX(3*n);
    std::getline( f, line );
    counter = 0;
    while( counter < 3*n &&  std::getline( f, line ) ){
        std::istringstream ss(line);
        for(auto i=0; i < numsPerLine; ++i){
            ss >> xi;
            pX(counter++) = xi;
        }
    }
    // Read x vector
    Vec xv(6*n);
    std::getline( f, line );
    counter = 0;
    while( counter < 6*n &&  std::getline( f, line ) ){
        std::istringstream ss(line);
        for(auto i=0; i < numsPerLine; ++i){
            ss >> xi;
            xv(counter++) = xi;
        }
    }

    f.close();
    SimulationState sims(n,ns,s,g,b,xv,pX,ip,l,en,rg);
    return sims;
}

void SimulationState::writeToFile(std::string file){
    // Overwrite any existing file
    ofstream f(file.c_str());
    f << "# !!! Auto-generated file DO NOT EDIT !!! #" << std::endl
      << "# Number of particles" << std::endl
      << N << std::endl
      << "# Latest output VTK file suffix" << std::endl
      << nameSuffix << std::endl
      << "# Number of time steps completed" << std::endl
      << step << std::endl
      << "# Gamma or the FvK number" << std::endl
      << gamma << std::endl
      << "# Beta or 1/Temperature" << std::endl
      << beta << std::endl
      << "# State of the random number engine" << std::endl
      << engine << std::endl
      << "# State of the random number generator" << std::endl
      << rng << std::endl
      << "# Initial nearest neighbors" << std::endl;

    size_t numLines = N/numsPerLine;
    size_t startIdx = 0;

    // Write nearest neighbors
    for(auto i = 0; i < numLines; ++i){
        startIdx = i*numsPerLine;
        for(auto j = 0; j < numsPerLine - 1; ++j){
            f << nn[ startIdx++ ] << " ";
        }
        f << nn[ startIdx++ ] << std::endl;
    }
    if( N % numsPerLine > 0){
        for( ; startIdx < N - 1; ){
            f << nn[ startIdx++ ] << " ";
        }
        f << nn[ startIdx ] << std::endl;
    }

    // Write initial position vector
    f << "# Particle positions at step 0" << std::endl;
    numLines = (3*N)/numsPerLine;
    for(auto i = 0; i < numLines; ++i){
        startIdx = i*numsPerLine;
        for(auto j = 0; j < numsPerLine - 1; ++j){
            f << initPos( startIdx++ ) << " ";
        }
        f << initPos( startIdx++ ) << std::endl;
    }
    if( 3*N % numsPerLine > 0){
        for( ; startIdx < 3*N - 1; ){
            f << initPos( startIdx++ ) << " ";
        }
        f << initPos( startIdx ) << std::endl;
    }

    // Write prevX vector
    f << "# Particle positions of previous step" << std::endl;
    numLines = (3*N)/numsPerLine;
    for(auto i = 0; i < numLines; ++i){
        startIdx = i*numsPerLine;
        for(auto j = 0; j < numsPerLine - 1; ++j){
            f << prevX( startIdx++ ) << " ";
        }
        f << prevX( startIdx++ ) << std::endl;
    }
    if( 3*N % numsPerLine > 0){
        for( ; startIdx < 3*N - 1; ){
            f << prevX( startIdx++ ) << " ";
        }
        f << prevX( startIdx ) << std::endl;
    }

    // Write x vector
    f << "# Particle positions and rotation vectors of latest step" << std::endl;
    numLines = (6*N)/numsPerLine;
    for(auto i = 0; i < numLines; ++i){
        startIdx = i*numsPerLine;
        for(auto j = 0; j < numsPerLine - 1; ++j){
            f << x( startIdx++ ) << " ";
        }
        f << x( startIdx++ ) << std::endl;
    }
    if( 6*N % numsPerLine > 0){
        for( ; startIdx < 6*N - 1; ){
            f << x( startIdx++ ) << " ";
        }
        f << x( startIdx ) << std::endl;
    }
    f.close();
}

}
