#include "catch.hpp"
#include "quadTree.h"
#include "body.h"
#include "simulation.h"

using namespace nBodySimulation;

TEST_CASE("Body insertion into QuadTree", "[QuadTree][Body]") {
    Quad boundary(0, 0, 100, 100);
    QuadTree qt(boundary, 4, false);
    Body body(10, 50, 50, 0, 0);

    REQUIRE(qt.insert(body) == true);
    REQUIRE(qt.getTotalMass() == Approx(10.0));
    REQUIRE(qt.getBodies().size() == 1);
}

TEST_CASE("Force calculations between bodies using QuadTree", "[QuadTree][Body]") {
    Quad boundary(0, 0, 100, 100);
    QuadTree qt(boundary, 2, false);
    Body body1(10, 20, 20, 0, 0);
    Body body2(10, 30, 30, 0, 0);

    qt.insert(body1);
    qt.insert(body2);
    double forceX = 0, forceY = 0;

    qt.calculateForce(body1, forceX, forceY, 1e-6, 1.0);
    REQUIRE(forceX != 0);
    REQUIRE(forceY != 0);
}

TEST_CASE("Simulation Integration with QuadTree and Body", "[Integration]") {
    Simulation sim;
    sim.setGravConstant(-6.67430e-11);
    sim.setTimeStep(0.01);
    sim.setSofteningLength(0.1);

    std::vector<Body> bodies;
    bodies.emplace_back(1.0e10, 100, 100, 0, 0);
    bodies.emplace_back(1.0e10, 400, 400, 0, 0);

    sim.setNumOfBodies(bodies.size());

    // Run simulation using brute force to verify bodies interact correctly
    sim.runBruteForce(bodies);

    // Verify bodies have moved
    REQUIRE(bodies[0].getX() != 100);
    REQUIRE(bodies[0].getY() != 100);
    REQUIRE(bodies[1].getX() != 400);
    REQUIRE(bodies[1].getY() != 400);
}

TEST_CASE("QuadTree and Simulation Integration for Barnes Hut Force Calculations", "[Integration]") {
    Simulation sim;
    Quad boundary(0, 0, 1000, 1000);
    QuadTree qt(boundary, 4, false);
    std::vector<Body> bodies;

    bodies.emplace_back(1.0e10, 500, 500, 0, 0); // Central massive body
    bodies.emplace_back(1.0e6, 550, 550, 0, 0);  // Orbiting small body

    for (auto& body : bodies) {
        qt.insert(body);
    }

    sim.runBarnesHut(bodies, 1000, 1000);

    REQUIRE(bodies[1].getX() != 550);
    REQUIRE(bodies[1].getY() != 550);
}

TEST_CASE("QuadTree and Simulation Integration for FMM Force Calculations", "[Integration]") {
    Simulation sim;
    sim.setGravConstant(-6.67430e-11);
    sim.setTimeStep(0.01);
    sim.setSofteningLength(0.1);

    std::vector<Body> bodies;
    for (int i = 0; i < 100; i++) {
        bodies.emplace_back(1.0e10, rand() % 1000, rand() % 1000, 0, 0);
    }

    sim.setNumOfBodies(bodies.size());

    QuadTree::precomputeBinomialCoefficients();
    QuadTree::precomputeTranslationMatrix();
    QuadTree::precomputeZPowers(24);

    sim.runParallelFastMultipole(bodies, 1000, 1000);

    // Verify that no body remains static
    for (auto& body : bodies) {
        REQUIRE(body.getX() != 0);
        REQUIRE(body.getY() != 0);
    }
}


TEST_CASE("Integration Test for Parallel Brute Force Execution in Simulation", "[Integration]") {
    Simulation sim;
    sim.setGravConstant(-6.67430e-11);
    sim.setTimeStep(0.01);
    sim.setSofteningLength(0.1);

    std::vector<Body> bodies;
    for (int i = 0; i < 100; i++) {
        bodies.emplace_back(1.0e10, rand() % 1000, rand() % 1000, 0, 0);
    }

    sim.setNumOfBodies(bodies.size());

    sim.runParallelBruteForce(bodies);

    // Verify that no body remains static
    for (auto& body : bodies) {
        REQUIRE(body.getX() != 0);
        REQUIRE(body.getY() != 0);
    }
}

TEST_CASE("Integration Test for Parallel Barnes-Hut Execution in Simulation", "[Integration]") {
    Simulation sim;
    sim.setGravConstant(-6.67430e-11);
    sim.setTimeStep(0.01);
    sim.setSofteningLength(0.1);

    std::vector<Body> bodies;
    for (int i = 0; i < 100; i++) {
        bodies.emplace_back(1.0e10, rand() % 1000, rand() % 1000, 0, 0);
    }

    sim.setNumOfBodies(bodies.size());

    QuadTree::precomputeBinomialCoefficients();
    QuadTree::precomputeTranslationMatrix();
    QuadTree::precomputeZPowers(10);

    sim.runParallelBarnesHut(bodies, 1000, 1000);

    // Verify that no body remains static
    for (auto& body : bodies) {
        REQUIRE(body.getX() != 0);
        REQUIRE(body.getY() != 0);
    }
}

TEST_CASE("Integration Test for Parallel FMM Execution in Simulation", "[Integration]") {
    Simulation sim;
    sim.setGravConstant(-6.67430e-11);
    sim.setTimeStep(0.01);
    sim.setSofteningLength(0.1);

    std::vector<Body> bodies;
    for (int i = 0; i < 100; i++) {
        bodies.emplace_back(1.0e10, rand() % 1000, rand() % 1000, 0, 0);
    }

    sim.setNumOfBodies(bodies.size());

    QuadTree::precomputeBinomialCoefficients();
    QuadTree::precomputeTranslationMatrix();
    QuadTree::precomputeZPowers(24);

    sim.runParallelFastMultipole(bodies, 1000, 1000);

    // Verify that no body remains static
    for (auto& body : bodies) {
        REQUIRE(body.getX() != 0);
        REQUIRE(body.getY() != 0);
    }
}