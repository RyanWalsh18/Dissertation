#include "quadTree.h"
#include "body.h"
#include <iostream>
#include <cmath>
#include <omp.h>

namespace nBodySimulation
{

    // ----------------------Multipole struct---------------------
    MultipoleExpansion::MultipoleExpansion(int numCoefficients)
    {
        this->mass = 0.0;
        this->coefficients.resize(numCoefficients, std::complex<double>(0.0, 0.0));
    }
    MultipoleExpansion::MultipoleExpansion()
    {
        this->mass = 0.0;
    }

    // ------------------------Quad class------------------------
    Quad::Quad(double topLeftX, double topLeftY, double bottomRightX, double bottomRightY)
    {
        this->topLeftX = topLeftX;
        this->topLeftY = topLeftY;
        this->bottomRightX = bottomRightX;
        this->bottomRightY = bottomRightY;
    }

    double Quad::getTopLeftX() const { return topLeftX; }
    double Quad::getTopLeftY() const { return topLeftY; }
    double Quad::getBottomRightX() const { return bottomRightX; }
    double Quad::getBottomRightY() const { return bottomRightY; }

    bool Quad::contains(const Body &body)
    {
        return body.getX() >= this->topLeftX && body.getX() <= this->bottomRightX && body.getY() >= this->topLeftY && body.getY() <= this->bottomRightY;
    }

    // --------------------QuadTree class------------------------
    std::vector<std::vector<double>> QuadTree::binomialCoefficients;
    std::vector<std::vector<std::complex<double>>> QuadTree::translation_matrix;
    std::vector<std::vector<std::complex<double>>> QuadTree::z_powers;

    void QuadTree::precomputeBinomialCoefficients()
    {
        binomialCoefficients.resize(25, std::vector<double>(25, 0.0));

        for (int n = 0; n < 25; n++)
        {
            binomialCoefficients.at(n).at(0) = 1.0;
            binomialCoefficients.at(n).at(n) = 1.0;

            for (int k = 1; k < n; k++)
            {
                binomialCoefficients.at(n).at(k) = binomialCoefficients.at(n - 1).at(k - 1) + binomialCoefficients.at(n - 1).at(k);
            }
        }
    }

    void QuadTree::precomputeTranslationMatrix()
    {
        translation_matrix.resize(24, std::vector<std::complex<double>>(24, 0.0));

        for (int j = 0; j < 24; j++)
        {
            for (int k = 0; k < j; k++)
            {
                translation_matrix.at(j).at(k) = binomialCoefficients.at(j).at(k);
            }
        }
    }

    void QuadTree::precomputeZPowers(int maxDegree)
    {
        z_powers.resize(maxDegree + 1, std::vector<std::complex<double>>(maxDegree + 1));
        for (int i = 0; i <= maxDegree; i++)
        {
            z_powers.at(i).at(0) = 1.0;
            for (int j = 1; j <= maxDegree; j++)
            {
                z_powers.at(i).at(j) = z_powers.at(i).at(j - 1) * std::complex<double>(std::cos(2 * 3.142 * i / maxDegree), std::sin(2 * 3.142 * i / maxDegree));
            }
        }
    }

    QuadTree::QuadTree(Quad boundary, int capacity, int newLevel, QuadTree *treeRoot)
    {
        this->boundary = std::make_unique<Quad>(boundary.getTopLeftX(), boundary.getTopLeftY(), boundary.getBottomRightX(), boundary.getBottomRightY());
        this->capacity = capacity;
        this->outer_coefficients = MultipoleExpansion(24);
        this->inner_coefficients = MultipoleExpansion(24);
        this->root = false;
        this->parent = nullptr;
        this->level = newLevel;
        this->treeRoot = treeRoot;
    }

    QuadTree::QuadTree(Quad boundary, int capacity, bool root)
    {
        this->boundary = std::make_unique<Quad>(boundary.getTopLeftX(), boundary.getTopLeftY(), boundary.getBottomRightX(), boundary.getBottomRightY());
        this->capacity = capacity;
        this->outer_coefficients = MultipoleExpansion(24);
        this->inner_coefficients = MultipoleExpansion(24);
        this->root = true;
        this->parent = nullptr;
        this->level = 0;
        this->treeRoot = this;
    }

    int QuadTree::getCapacity() const { return capacity; }
    bool QuadTree::isDivided() const { return divided; }
    double QuadTree::getTotalMass() const { return totalMass; }
    double QuadTree::getComX() const { return comX; }
    double QuadTree::getComY() const { return comY; }
    bool QuadTree::isRoot() const { return root; }
    int QuadTree::getLevel() const { return level; }
    double QuadTree::getPotential() const { return potential; }
    const std::vector<Body>& QuadTree::getBodies() const { return bodies; }

    QuadTree* QuadTree::getNorthwest() const { return northwest.get(); }
    QuadTree* QuadTree::getNortheast() const { return northeast.get(); }
    QuadTree* QuadTree::getSouthwest() const { return southwest.get(); }
    QuadTree* QuadTree::getSoutheast() const { return southeast.get(); }
    
    void QuadTree::subdivide()
    {
        double xMid = (boundary->getTopLeftX() + boundary->getBottomRightX()) / 2.0;
        double yMid = (boundary->getTopLeftY() + boundary->getBottomRightY()) / 2.0;

        this->northwest = std::make_unique<QuadTree>(Quad(boundary->getTopLeftX(), boundary->getTopLeftY(), xMid, yMid), capacity, this->level + 1, this->treeRoot);
        this->northwest->parent = this;

        this->northeast = std::make_unique<QuadTree>(Quad(xMid, boundary->getTopLeftY(), boundary->getBottomRightX(), yMid), capacity, this->level + 1, this->treeRoot);
        this->northeast->parent = this;

        this->southwest = std::make_unique<QuadTree>(Quad(boundary->getTopLeftX(), yMid, xMid, boundary->getBottomRightY()), capacity, this->level + 1, this->treeRoot);
        this->southwest->parent = this;

        this->southeast = std::make_unique<QuadTree>(Quad(xMid, yMid, boundary->getBottomRightX(), boundary->getBottomRightY()), capacity, this->level + 1, this->treeRoot);
        this->southeast->parent = this;

        this->divided = true;
    }

    void QuadTree::updateMassAndCOM(const Body &body)
    {
        double newTotalMass = this->totalMass + body.getMass();

        this->comX = (this->comX * this->totalMass + body.getX() * body.getMass()) / newTotalMass;
        this->comY = (this->comY * this->totalMass + body.getY() * body.getMass()) / newTotalMass;

        this->totalMass = newTotalMass;
    }

    bool QuadTree::insert(const Body &body)
    {
        if (!(this->boundary->contains(body)))
        {
            return false;
        }

        updateMassAndCOM(body);

        int bodiesInNode = this->bodies.size();

        //if current node has space/is empty and is a leaf node -> insert body
        if (bodiesInNode < this->capacity && !this->divided)
        {
            bodies.push_back(body);
            return true;
        }
        else
        {
            // currnet node is a leaf node and is full -> subdivide current node and reinsert current body and 'inserting' body further down tree
            if (!this->divided)
            {
                subdivide();

                std::vector<Body> existingBodies = std::move(bodies); // Move bodies to a temporary container
                bodies.clear();                                       // Clear the current node's bodies as they will be redistributed

                for (const Body &existingBody : existingBodies)
                {
                    if (!northwest->insert(existingBody) && !northeast->insert(existingBody) &&
                        !southwest->insert(existingBody) && !southeast->insert(existingBody))
                    {
                        std::cout << "Failed to insert existing body into any child node." << std::endl;
                        return false;
                    }
                }
            }
            // current is not a leaf node -> recusivley pass to children
            if (this->northwest->insert(body) || this->northeast->insert(body) || this->southwest->insert(body) || this->southeast->insert(body))
            {
                return true;
            }
            else{
                return false;
            }
        }

        return false;
    }

    void QuadTree::calculateForce(const Body &body, double &forceX, double &forceY, double softeningLen, double gravConst)
    {
        // if current is an external node -> caluclate the force directly
        if (this->divided == false)
        {

            for (Body &currentBody : this->bodies)
            {
                if (&currentBody != &body)
                {
                    double dx = currentBody.getX() - body.getX();
                    double dy = currentBody.getY() - body.getY();
                    double dist = sqrt(dx * dx + dy * dy + softeningLen * softeningLen);
                    double forceMagnitude = gravConst * currentBody.getMass() / (dist * dist * dist);
                    forceX += forceMagnitude * dx;
                    forceY += forceMagnitude * dy;
                }
            }
        }
        else
        {
            // calculate the s/d ratio. if ratio < theta treat node as single body using com values
            double dx = this->comX - body.getX();
            double dy = this->comY - body.getY();
            double dist = sqrt(dx * dx + dy * dy);

            if ((this->boundary->getBottomRightX() - this->boundary->getTopLeftX()) / dist < 0.5)
            {
                double forceMagnitude = gravConst * this->totalMass / (dist * dist + softeningLen * softeningLen);
                forceX += forceMagnitude * dx / dist;
                forceY += forceMagnitude * dy / dist;
            }
            else
            {
                // run recursively on the children
                if (this->northwest)
                {
                    this->northwest->calculateForce(body, forceX, forceY, softeningLen, gravConst);
                }
                if (northeast)
                {
                    this->northeast->calculateForce(body, forceX, forceY, softeningLen, gravConst);
                }
                if (southwest)
                {
                    this->southwest->calculateForce(body, forceX, forceY, softeningLen, gravConst);
                }
                if (southeast)
                {
                    this->southeast->calculateForce(body, forceX, forceY, softeningLen, gravConst);
                }
            }
        }
    }

    void QuadTree::compute_leaf_outer_expansion()
    {
        double softeningLen = 1e-6;

        double x_center = (this->boundary->getTopLeftX() + this->boundary->getBottomRightX()) / 2;
        double y_center = (this->boundary->getTopLeftY() + this->boundary->getBottomRightY()) / 2;

        std::complex<double> z_center(x_center, y_center);

        MultipoleExpansion &outer = this->outer_coefficients;

        double total_mass = 0.0;
        int i = 0;
        for (i = 0; i < this->bodies.size(); i++)
        {
            std::complex<double> z(this->bodies.at(i).getX(), this->bodies.at(i).getY());

            std::complex<double> z_relative = z - z_center;
            double z_relative_abs = std::sqrt(std::norm(z_relative) + softeningLen * softeningLen);

            double mass = this->bodies.at(i).getMass();

            total_mass += mass;

            for (int e = 1; e <= 24; e++)
            {
                
        
                outer.coefficients.at(e - 1) += mass * z_powers.at(i).at(e) / std::pow(z_relative_abs, 2 * e);
                
            }
        }

        outer.mass = total_mass;
    }

    void QuadTree::outer_shift(MultipoleExpansion &parent_expansion, const MultipoleExpansion &child_expansion, std::complex<double> &z_diff)
    {

        parent_expansion.mass += child_expansion.mass;

        int j=0;
        for (j = 0; j < 24; j++)
        {
            std::complex<double> sum = 0.0;
            int k = 0;
            for (k = 0; k < 24; k++)
            {
                sum += translation_matrix.at(j).at(k) * child_expansion.coefficients.at(k) * std::pow(z_diff, j - k);
            }
            
            
            parent_expansion.coefficients.at(j) += sum;
            
        }
    }

    void QuadTree::build_outer_expansions()
    {
        // if node is a leaf -> compute outer expansion directly
        if (!(this->divided))
        {
            compute_leaf_outer_expansion();
        }
        else
        {
            // post order traversal
            this->northeast->build_outer_expansions();
            this->northwest->build_outer_expansions();
            this->southeast->build_outer_expansions();
            this->southwest->build_outer_expansions();

            MultipoleExpansion aggregate_expansion(24);

            if (this->northeast)
            {
                std::complex<double> z_diff = std::complex<double>(this->northeast->comX, this->northeast->comY) - std::complex<double>(this->comX, this->comY);
                outer_shift(aggregate_expansion, this->northeast->outer_coefficients, z_diff);
            }
            if (this->northwest)
            {
                std::complex<double> z_diff = std::complex<double>(this->northwest->comX, this->northwest->comY) - std::complex<double>(this->comX, this->comY);
                outer_shift(aggregate_expansion, this->northwest->outer_coefficients, z_diff);
            }
            if (this->southeast)
            {
                std::complex<double> z_diff = std::complex<double>(this->southeast->comX, this->southeast->comY) - std::complex<double>(this->comX, this->comY);
                outer_shift(aggregate_expansion, this->southeast->outer_coefficients, z_diff);
            }
            if (this->southwest)
            {
                std::complex<double> z_diff = std::complex<double>(this->southwest->comX, this->southwest->comY) - std::complex<double>(this->comX, this->comY);
                outer_shift(aggregate_expansion, this->southwest->outer_coefficients, z_diff);
            }

            this->outer_coefficients = aggregate_expansion;
        }
    }

    void QuadTree::inner_shift(MultipoleExpansion &inner_expansion, const std::complex<double> &z_diff)
    {
        std::vector<std::complex<double> > shifted_coefficients(inner_expansion.coefficients.size());

        for (int n = 0; n < inner_expansion.coefficients.size(); n++)
        {
            shifted_coefficients.at(n) = inner_expansion.coefficients.at(n);

            for (int k = 1; k <= n; k++)
            {
                shifted_coefficients.at(n) += translation_matrix.at(n).at(k) * inner_expansion.coefficients.at(n - k) * std::pow(z_diff, k);
            }
        }

        inner_expansion.coefficients = shifted_coefficients;
    }

    void QuadTree::get_neighbors_at_level(std::vector<QuadTree *> &neighbors, int target_level)
    {
        if (this->level == target_level)
        {
            neighbors.push_back(this);
        }

        else if (this->level < target_level && this->divided)
        {
            this->northwest->get_neighbors_at_level(neighbors, target_level);
            this->northeast->get_neighbors_at_level(neighbors, target_level);
            this->southwest->get_neighbors_at_level(neighbors, target_level);
            this->southeast->get_neighbors_at_level(neighbors, target_level);
        }
    }

    bool QuadTree::is_direct_neighbor(QuadTree *node)
    {
        return (node == this->northwest.get() || node == this->northeast.get() ||
                node == this->southwest.get() || node == this->southeast.get());
    }

    void QuadTree::compute_interaction_set()
    {
        this->interaction_set.clear();

        std::vector<QuadTree *> parent_neighbors;
        this->treeRoot->get_neighbors_at_level(parent_neighbors, this->level - 1);

        for (QuadTree *neighbor : parent_neighbors)
        {
            if (neighbor->divided)
            {
                if (!is_direct_neighbor(neighbor->northwest.get()))
                {
                    interaction_set.push_back(neighbor->northwest.get());
                }
                if (!is_direct_neighbor(neighbor->northeast.get()))
                {
                    interaction_set.push_back(neighbor->northeast.get());
                }
                if (!is_direct_neighbor(neighbor->southwest.get()))
                {
                    interaction_set.push_back(neighbor->southwest.get());
                }
                if (!is_direct_neighbor(neighbor->southeast.get()))
                {
                    interaction_set.push_back(neighbor->southeast.get());
                }
            }
        }
    }

    void QuadTree::outer_to_inner(MultipoleExpansion &inner_expansion, const MultipoleExpansion &outer_expansion, const std::complex<double> &z_diff)
    {
        inner_expansion.mass = outer_expansion.mass;

        for (int n = 0; n < outer_expansion.coefficients.size(); n++)
        {
            inner_expansion.coefficients.at(n) = std::complex<double>(0.0, 0.0);

            for (int k = 0; k <= n; k++)
            {
                inner_expansion.coefficients.at(n) += binomialCoefficients.at(n).at(k) * outer_expansion.coefficients.at(k) * std::pow(-z_diff, n - k);
            }

            inner_expansion.coefficients.at(n) *= std::pow(-1, n);
        }
    }

    void QuadTree::add_expansions(MultipoleExpansion &target_expansion, const MultipoleExpansion &source_expansion)
    {
        target_expansion.mass += source_expansion.mass;

        for (int n = 0; n < target_expansion.coefficients.size(); n++)
        {
            target_expansion.coefficients.at(n) += source_expansion.coefficients.at(n);
        }
    }

    void QuadTree::get_direct_neighbors(std::vector<QuadTree *> &neighbors)
    {
        QuadTree *current = this;
        QuadTree *parent = current->parent;

        while (parent != nullptr)
        {
            double currentX = current->boundary->getTopLeftX() + (current->boundary->getBottomRightX() - current->boundary->getTopLeftX()) / 2;
            double currentY = current->boundary->getTopLeftY() + (current->boundary->getBottomRightY() - current->boundary->getTopLeftY()) / 2;

            double parentX = parent->boundary->getTopLeftX() + (parent->boundary->getBottomRightX() - parent->boundary->getTopLeftX()) / 2;
            double parentY = parent->boundary->getTopLeftY() + (parent->boundary->getBottomRightY() - parent->boundary->getTopLeftY()) / 2;

            if (currentX < parentX && currentY < parentY)
            {
                // Current node is the southwest child
                if (parent->northwest && parent->northwest->divided)
                {
                    neighbors.push_back(parent->northwest->southwest.get());
                    neighbors.push_back(parent->northwest->southeast.get());
                }
                if (parent->northeast && parent->northeast->divided)
                {
                    neighbors.push_back(parent->northeast->southwest.get());
                }
                if (parent->southeast && parent->southeast->divided)
                {
                    neighbors.push_back(parent->southeast->northwest.get());
                }
            }
            else if (currentX < parentX && currentY >= parentY)
            {
                // Current node is the northwest child
                if (parent->southwest && parent->southwest->divided)
                {
                    neighbors.push_back(parent->southwest->northeast.get());
                }
                if (parent->southeast && parent->southeast->divided)
                {
                    neighbors.push_back(parent->southeast->northwest.get());
                    neighbors.push_back(parent->southeast->northeast.get());
                }
                if (parent->northeast && parent->northeast->divided)
                {
                    neighbors.push_back(parent->northeast->southwest.get());
                }
            }
            else if (currentX >= parentX && currentY < parentY)
            {
                // Current node is the southeast child
                if (parent->southwest && parent->southwest->divided)
                {
                    neighbors.push_back(parent->southwest->northeast.get());
                }
                if (parent->northwest && parent->northwest->divided)
                {
                    neighbors.push_back(parent->northwest->southeast.get());
                    neighbors.push_back(parent->northwest->northeast.get());
                }
                if (parent->northeast && parent->northeast->divided)
                {
                    neighbors.push_back(parent->northeast->southwest.get());
                }
            }
            else
            {
                // Current node is the northeast child
                if (parent->northwest && parent->northwest->divided)
                {
                    neighbors.push_back(parent->northwest->southeast.get());
                }
                if (parent->southwest && parent->southwest->divided)
                {
                    neighbors.push_back(parent->southwest->northeast.get());
                }
                if (parent->southeast && parent->southeast->divided)
                {
                    neighbors.push_back(parent->southeast->northwest.get());
                    neighbors.push_back(parent->southeast->southwest.get());
                }
            }

            current = parent;
            parent = current->parent;
        }
    }

    void QuadTree::build_inner_expansions()
    {
        // if this is the root -> initialise parent inner_expansion to zero
        if (this->root)
        {
            this->inner_coefficients = MultipoleExpansion(24);
        }
        else
        {
            // shift the parents inner_coefficients to the child
            this->inner_coefficients = parent->inner_coefficients;
            std::complex<double> z_diff = std::complex<double>(this->comX, this->comY) - std::complex<double>(this->parent->comX, this->parent->comY);

            inner_shift(this->inner_coefficients, z_diff);

            // for each node in the interaction set, convert the outer coefficients into the inner coefficients
            compute_interaction_set();

            int i = 0;
            for (i = 0; i < interaction_set.size(); i++)
            {
                QuadTree *node = interaction_set.at(i);

                MultipoleExpansion temp_inner(24);
                std::complex<double> z_diff = std::complex<double>(this->comX, this->comY) - std::complex<double>(node->comX, node->comY);
                outer_to_inner(temp_inner, node->outer_coefficients, z_diff);

        
                
                add_expansions(this->inner_coefficients, temp_inner);
                
                
            }
        }

        if (this->divided)
        {
            if (this->northeast)
            {
                this->northeast->build_inner_expansions();
            }
            if (this->northwest)
            {
                this->northwest->build_inner_expansions();
            }
            if (this->southeast)
            {
                this->southeast->build_inner_expansions();
            }
            if (this->southwest)
            {
                this->southwest->build_inner_expansions();
            }
        }
        else
        {
            double softeningLen = 1e-6;

            // Leaf node: directly compute inner coefficients from 8 direct neighbors
            std::vector<QuadTree *> neighbors;
            get_direct_neighbors(neighbors);

            int i = 0;
            for (i = 0; i < neighbors.size(); i++)
            {
                QuadTree *neighbor = neighbors.at(i);
                if (!neighbor->divided)
                {
                    for (const Body &body : neighbor->bodies)
                    {
                        std::complex<double> z_diff = std::complex<double>(body.getX(), body.getY()) - std::complex<double>(this->comX, this->comY);
                        double dist = std::sqrt(std::norm(z_diff) + softeningLen * softeningLen);
                        double potential = std::log(dist);

             
                        
                        this->inner_coefficients.mass += body.getMass();
                        this->inner_coefficients.coefficients.at(0) += body.getMass() * potential;
                        

                        std::complex<double> z_pow = 1.0;
                        for (int n = 1; n < this->inner_coefficients.coefficients.size(); n++)
                        {
                            z_pow *= z_diff;
           
                            
                            this->inner_coefficients.coefficients.at(n) += body.getMass() * z_pow / std::pow(dist, 2 * n);
                            
                        }
                    }
                }
            }
        }

        if (!this->root)
        {
            std::complex<double> z_diff = std::complex<double>(this->comX, this->comY) - std::complex<double>(this->parent->comX, this->parent->comY);
            this->potential = this->inner_coefficients.coefficients.at(0).real();

            std::complex<double> z_pow = 1.0;
            for (int n = 1; n < this->inner_coefficients.coefficients.size(); n++)
            {
                z_pow *= z_diff;
                this->potential += (this->inner_coefficients.coefficients.at(n) * z_pow).real();
            }
        }
        else
        {
            this->potential = this->inner_coefficients.coefficients.at(0).real();
        }
    }


    void QuadTree::calculateForcesFMM(Body &body, double &forceX, double &forceY, double softeningLen, double gravConst) {
        const double epsilon = 1e-10;  // A small constant to prevent division by zero

        if (!this->divided) {
            // Direct force calculation for leaf nodes
            for (const Body &otherBody : this->bodies) {
                if (&body != &otherBody) {
                    double dx = otherBody.getX() - body.getX();
                    double dy = otherBody.getY() - body.getY();
                    double distSq = dx * dx + dy * dy + softeningLen * softeningLen;
                    double dist = std::sqrt(distSq);

                    // Ensuring that we do not divide by zero
                    if (distSq > epsilon) {
                        double forceMagnitude = gravConst * body.getMass() * otherBody.getMass() / distSq;
                        forceX += forceMagnitude * dx / dist;
                        forceY += forceMagnitude * dy / dist;
                    }
                }
            }
        } else {
            // Use multipole expansion for internal nodes
            double dx = this->comX - body.getX();
            double dy = this->comY - body.getY();
            double distSq = dx * dx + dy * dy;
            double dist = std::sqrt(distSq + softeningLen * softeningLen);
            double sideLength = this->boundary->getBottomRightX() - this->boundary->getTopLeftX();

            if ((sideLength / dist) < 0.1 && distSq > epsilon) {
                std::complex<double> z_diff(dx, dy);
                std::complex<double> z_pow(1.0, 0.0);

                for (int n = 1; n < this->inner_coefficients.coefficients.size(); ++n) {
                    z_pow *= z_diff; // Increment the power of z_diff
                    std::complex<double> expansionTerm = this->inner_coefficients.coefficients[n] * z_pow;
                    double forceTerm = gravConst * body.getMass() * (n / std::pow(dist, n + 1));
                    forceX += expansionTerm.real() * forceTerm;
                    forceY += expansionTerm.imag() * forceTerm;
                }
            } else {
                // Recurse on children for closer nodes
                if (northwest) northwest->calculateForcesFMM(body, forceX, forceY, softeningLen, gravConst);
                if (northeast) northeast->calculateForcesFMM(body, forceX, forceY, softeningLen, gravConst);
                if (southwest) southwest->calculateForcesFMM(body, forceX, forceY, softeningLen, gravConst);
                if (southeast) southeast->calculateForcesFMM(body, forceX, forceY, softeningLen, gravConst);
            }
        }
    }


}