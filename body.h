#ifndef BODY_H
#define BODY_H

namespace nBodySimulation
{

    class Body
    {
    private:
        double mass;
        double x, y;
        double vx, vy;
        double ax, ay;

    public:
        // Constructor
        Body();

        Body(double bodyMass, double xPos, double yPos, double xVel, double yVel, double xAcc, double yAcc);

        Body(double bodyMass, double xPos, double yPos, double xVel, double yVel);

        friend bool operator==(const Body &lhs, const Body &rhs);

        friend bool operator!=(const Body &lhs, const Body &rhs);

        // Getter functions
        double getMass() const;
        double getX() const;
        double getY() const;
        double getVX() const;
        double getVY() const;
        double getAX() const;
        double getAY() const;

        // Setter functions
        void setMass(double newMass);
        void setX(double xNew);
        void setY(double yNew);
        void setVX(double vxNew);
        void setVY(double vyNew);
        void setAX(double axNew);
        void setAY(double ayNew);

        void setPosition(double xNew, double yNew);
        void setVelocity(double vxNew, double vyNew);
        void setAcceleration(double axNew, double ayNew);
    };
}

#endif
