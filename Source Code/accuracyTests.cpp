#include "body.h"
#include "quadTree.h"
#include "simulation.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "accuracyTests.h"

using namespace std;
using namespace nBodySimulation;

void brute_force(double time_period, int numOfOrbits, double dt)
{
    std::vector<nBodySimulation::Body> bodies;

    Body body1 = Body(10000, 160, 120, 0, 0);
    Body body2 = Body(1, 260, 120, 0, 10.0);

    bodies.push_back(body1);
    bodies.push_back(body2);

    nBodySimulation::Simulation simulation(0.0, time_period * numOfOrbits, dt, 1e-6, -(1));

    int finished = 0;

    ofstream brute_force_file("brute_force_orbit.csv");
    brute_force_file << "Time,x,y,vx,vy\n";

    while (finished == 0)
    {
        if (simulation.getCurrentSimTime() < simulation.getSimEndTime())
        {
            brute_force_file << simulation.getCurrentSimTime() << "," << bodies.at(1).getX() << "," << bodies.at(1).getY() << "," << bodies.at(1).getVX() << "," << bodies.at(1).getVY() << "\n";
            simulation.runBruteForce(bodies);
        }
        else
        {
            brute_force_file.close();
            finished++;
        }
    }
}

void barnes_hut(double time_period, int numOfOrbits, double dt)
{
    std::vector<nBodySimulation::Body> bodies;

    Body body1 = Body(10000, 160, 120, 0, 0);
    Body body2 = Body(1, 260, 120, 0, 10.0);

    bodies.push_back(body1);
    bodies.push_back(body2);

    nBodySimulation::Simulation simulation(0.0, time_period * numOfOrbits, dt, 1e-6, -(1));

    int finished = 0;

    ofstream barnes_hut_file("barnes_hut_orbit.csv");
    barnes_hut_file << "Time,x,y,vx,vy\n";

    while (finished == 0)
    {
        if (simulation.getCurrentSimTime() < simulation.getSimEndTime())
        {
            barnes_hut_file << simulation.getCurrentSimTime() << "," << bodies.at(1).getX() << "," << bodies.at(1).getY() << "," << bodies.at(1).getVX() << "," << bodies.at(1).getVY() << "\n";
            simulation.runBarnesHut(bodies,300,300);
        }
        else
        {
            barnes_hut_file.close();
            finished++;
        }
    }
}

void fast_multipole(double time_period, int numOfOrbits, double dt){
    std::vector<nBodySimulation::Body> bodies;

    Body body1 = Body(10000, 160, 120, 0, 0);
    Body body2 = Body(1, 260, 120, 0, 10.0);

    bodies.push_back(body1);
    bodies.push_back(body2);

    nBodySimulation::Simulation simulation(0.0, time_period * numOfOrbits, dt, 1e-6, -(1));

    int finished = 0;

    ofstream fmm_file("fmm_orbit.csv");
    fmm_file << "Time,x,y,vx,vy\n";

    while (finished == 0)
    {
        QuadTree::precomputeBinomialCoefficients();
        QuadTree::precomputeTranslationMatrix();
        QuadTree::precomputeZPowers(24);

        if (simulation.getCurrentSimTime() < simulation.getSimEndTime())
        {
            fmm_file << simulation.getCurrentSimTime() << "," << bodies.at(1).getX() << "," << bodies.at(1).getY() << "," << bodies.at(1).getVX() << "," << bodies.at(1).getVY() << "\n";
            simulation.runFastMultipole(bodies,300,300);
        }
        else
        {
            fmm_file.close();
            finished++;
        }
    }
}

void parallel_brute_force(double time_period, int numOfOrbits, double dt){
    std::vector<nBodySimulation::Body> bodies;

    Body body1 = Body(10000, 160, 120, 0, 0);
    Body body2 = Body(1, 260, 120, 0, 10.0);

    bodies.push_back(body1);
    bodies.push_back(body2);

    nBodySimulation::Simulation simulation(0.0, time_period * numOfOrbits, dt, 1e-6, -(1));

    int finished = 0;

    ofstream parallel_brute_force_file("parallel_brute_force_orbit.csv");
    parallel_brute_force_file << "Time,x,y,vx,vy\n";

    while (finished == 0)
    {
        if (simulation.getCurrentSimTime() < simulation.getSimEndTime())
        {
            parallel_brute_force_file << simulation.getCurrentSimTime() << "," << bodies.at(1).getX() << "," << bodies.at(1).getY() << "," << bodies.at(1).getVX() << "," << bodies.at(1).getVY() << "\n";
            simulation.runParallelBruteForce(bodies);
        }
        else
        {
            parallel_brute_force_file.close();
            finished++;
        }
    }
}

void parallel_barnes_hut(double time_period, int numOfOrbits, double dt){
    std::vector<nBodySimulation::Body> bodies;

    Body body1 = Body(10000, 160, 120, 0, 0);
    Body body2 = Body(1, 260, 120, 0, 10.0);

    bodies.push_back(body1);
    bodies.push_back(body2);

    nBodySimulation::Simulation simulation(0.0, time_period * numOfOrbits, dt, 1e-6, -(1));

    int finished = 0;

    ofstream parallel_barnes_hut_file("parallel_barnes_hut_orbit.csv");
    parallel_barnes_hut_file << "Time,x,y,vx,vy\n";

    while (finished == 0)
    {
        if (simulation.getCurrentSimTime() < simulation.getSimEndTime())
        {
            parallel_barnes_hut_file << simulation.getCurrentSimTime() << "," << bodies.at(1).getX() << "," << bodies.at(1).getY() << "," << bodies.at(1).getVX() << "," << bodies.at(1).getVY() << "\n";
            simulation.runParallelBarnesHut(bodies,300,300);
        }
        else
        {
            parallel_barnes_hut_file.close();
            finished++;
        }
    }
}

void parallel_fast_multipole(double time_period, int numOfOrbits, double dt){
    std::vector<nBodySimulation::Body> bodies;

    Body body1 = Body(10000, 160, 120, 0, 0);
    Body body2 = Body(1, 260, 120, 0, 10.0);

    bodies.push_back(body1);
    bodies.push_back(body2);

    nBodySimulation::Simulation simulation(0.0, time_period * numOfOrbits, dt, 1e-6, -(1));

    int finished = 0;

    ofstream parallel_fmm_file("parallel_fmm_orbit.csv");
    parallel_fmm_file << "Time,x,y,vx,vy\n";

    while (finished == 0)
    {
        QuadTree::precomputeBinomialCoefficients();
        QuadTree::precomputeTranslationMatrix();
        QuadTree::precomputeZPowers(24);

        if (simulation.getCurrentSimTime() < simulation.getSimEndTime())
        {
            parallel_fmm_file << simulation.getCurrentSimTime() << "," << bodies.at(1).getX() << "," << bodies.at(1).getY() << "," << bodies.at(1).getVX() << "," << bodies.at(1).getVY() << "\n";
            simulation.runParallelFastMultipole(bodies,300,300);
        }
        else
        {
            parallel_fmm_file.close();
            finished++;
        }
    }
}

int main(int argc, char **argv){
    // M = 10000
    // r = 100
    // G = 1
    // Orbital time period (time for one orbit)= 60
    // dt = 0.01
    // orbital_velocity = 10.471

     // Constants
    double dt = 0.1;
    double time_period = 62.831;
    double ang_freq = 2 * M_PI / time_period;
    
    // Center of the orbit
    double centerX = 160;
    double centerY = 120;

    // Initial conditions of the orbiting body
    double initialX = 260;
    double initialY = 120;

    // Calculate the radius using the distance formula
    double r = sqrt(pow(initialX - centerX, 2) + pow(initialY - centerY, 2));
    
    // Calculate initial angle
    double initialTheta = atan2(initialY - centerY, initialX - centerX);

    // Initial position and velocity
    double t = 0;
    double x = centerX + r * cos(initialTheta);
    double y = centerY + r * sin(initialTheta);
    double vx = -r * ang_freq * sin(initialTheta);
    double vy = r * ang_freq * cos(initialTheta);

    int numOfOrbits = 5;

    ofstream benchmark_file("benchmark_orbit.csv");
    benchmark_file << "Time,x,y,vx,vy\n";

    while (t <= time_period*numOfOrbits) {
        benchmark_file << t << "," << x << "," << y << "," << vx << "," << vy << "\n";

        // Update time
        t += dt;

        // Calculate new angle
        double theta = initialTheta + ang_freq * t;

        // Update position
        x = centerX + r * cos(theta);
        y = centerY + r * sin(theta);

        // Update velocity
        vx = -r * ang_freq * sin(theta);
        vy = r * ang_freq * cos(theta);
    }

    benchmark_file.close();

    brute_force(time_period, numOfOrbits, dt);
    barnes_hut(time_period, numOfOrbits, dt);
    fast_multipole(time_period, numOfOrbits, dt);

    parallel_brute_force(time_period, numOfOrbits, dt);
    parallel_barnes_hut(time_period, numOfOrbits, dt);
    parallel_fast_multipole(time_period, numOfOrbits, dt);
}
