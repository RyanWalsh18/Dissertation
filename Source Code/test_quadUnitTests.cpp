#include "catch.hpp"
#include "quadTree.h"
#include "body.h"

using namespace nBodySimulation;

TEST_CASE("Quad initialises and reports dimensions correctly", "[Quad]") {
    Quad quad(1.0, 1.0, 5.0, 5.0);

    REQUIRE(quad.getTopLeftX() == 1.0);
    REQUIRE(quad.getTopLeftY() == 1.0);
    REQUIRE(quad.getBottomRightX() == 5.0);
    REQUIRE(quad.getBottomRightY() == 5.0);
}

TEST_CASE("Quad contains method works correctly", "[Quad]") {
    Quad quad(0.0, 0.0, 10.0, 10.0);
    Body insideBody(1.0, 5.0, 5.0, 0.0, 0.0);
    Body onEdgeBody(1.0, 10.0, 10.0, 0.0, 0.0);
    Body outsideBody(1.0, 15.0, 15.0, 0.0, 0.0);

    SECTION("Body inside the Quad") {
        REQUIRE(quad.contains(insideBody) == true);
    }

    SECTION("Body on the edge of the Quad") {
        REQUIRE(quad.contains(onEdgeBody) == true);
    }

    SECTION("Body outside the Quad") {
        REQUIRE(quad.contains(outsideBody) == false);
    }
}

TEST_CASE("Quad boundary conditions are handled correctly", "[Quad]") {
    Quad quad(0.0, 0.0, 10.0, 10.0);
    Body boundaryBody1(1.0, 0.0, 0.0, 0.0, 0.0); // Corner case
    Body boundaryBody2(1.0, 10.0, 10.0, 0.0, 0.0); // Opposite corner

    REQUIRE(quad.contains(boundaryBody1) == true);
    REQUIRE(quad.contains(boundaryBody2) == true);
}