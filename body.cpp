#include "body.h"
namespace nBodySimulation
{

    Body::Body()
    {
        this->mass = 0.0;
        this->x = 0.0;
        this->y = 0.0;
        this->vx = 0.0;
        this->vy = 0.0;
        this->ax = 0.0;
        this->ay = 0.0;
    }
    // Constructor without using an initializer list
    Body::Body(double bodyMass, double xPos, double yPos, double xVel, double yVel, double xAcc, double yAcc)
    {
        // Assigning values to member variables in the constructor body
        mass = bodyMass;
        x = xPos;
        y = yPos;
        vx = xVel;
        vy = yVel;
        ax = xAcc;
        ay = yAcc;
    }

    // Constructor without acceleration parameters (defaults to 0)
    Body::Body(double bodyMass, double xPos, double yPos, double xVel, double yVel)
    {
        mass = bodyMass;
        x = xPos;
        y = yPos;
        vx = xVel;
        vy = yVel;
        ax = 0;
        ay = 0;
    }

    bool operator==(const Body &lhs, const Body &rhs)
    {
        return lhs.mass == rhs.mass &&
               lhs.x == rhs.x &&
               lhs.y == rhs.y &&
               lhs.vx == rhs.vx &&
               lhs.vy == rhs.vy &&
               lhs.ax == rhs.ax &&
               lhs.ay == rhs.ay;
    }

    bool operator!=(const Body &lhs, const Body &rhs)
    {
        return lhs.mass != rhs.mass ||
               lhs.x != rhs.x ||
               lhs.y != rhs.y ||
               lhs.vx != rhs.vx ||
               lhs.vy != rhs.vy ||
               lhs.ax != rhs.ax ||
               lhs.ay != rhs.ay;
    }

    // Getter functions
    double Body::getMass() const { return mass; }
    double Body::getX() const { return x; }
    double Body::getY() const { return y; }
    double Body::getVX() const { return vx; }
    double Body::getVY() const { return vy; }
    double Body::getAX() const { return ax; }
    double Body::getAY() const { return ay; }

    // Setter functions
    void Body::setMass(double bodyMass) { mass = bodyMass; }
    void Body::setX(double xPos) { x = xPos; }
    void Body::setY(double yPos) { y = yPos; }
    void Body::setVX(double xVel) { vx = xVel; }
    void Body::setVY(double yVel) { vy = yVel; }
    void Body::setAX(double xAcc) { ax = xAcc; }
    void Body::setAY(double yAcc) { ay = yAcc; }

    void Body::setPosition(double xPos, double yPos)
    {
        x = xPos;
        y = yPos;
    }

    void Body::setVelocity(double xVel, double yVel)
    {
        vx = xVel;
        vy = yVel;
    }

    void Body::setAcceleration(double xAcc, double yAcc)
    {
        ax = xAcc;
        ay = yAcc;
    }
}
