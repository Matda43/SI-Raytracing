#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>

using namespace cv;
using namespace std;

struct Vec3
{
    float x, y, z;

    void Debug() const
    {
        std::cout << "(" << x << "," << y << "," << z << ")";
    }

    Vec3 operator*(const float f) const
    {
        return Vec3{ x * f, y * f, z * f };
    }

    Vec3 operator/(const float f) const
    {
        return Vec3{ x / f, y / f, z / f };
    }

    Vec3 operator*(const Vec3 v) const
    {
        return Vec3{ x * v.x, y * v.y, z * v.z };
    }

    Vec3 operator+(const Vec3 v) const
    {
        return Vec3{ x + v.x, y + v.y, z + v.z };
    }

    Vec3 operator-(const Vec3 v) const
    {
        return Vec3{ x - v.x, y - v.y, z - v.z };
    }

    float normSquared() const
    {
        return x * x + y * y + z * z;
    }

    Vec3 unitVector() const
    {
        const float norm = sqrt(normSquared());
        return Vec3{ x / norm, y / norm, z / norm };
    }

    float scalarProduct(Vec3 v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }

    Vec3 getNormal(Vec3 v) {
        Vec3 n = Vec3{ x - v.x, y - v.y, z - v.z };
        float n2 = n.normSquared();
        return n / sqrt(n2);
    }
};

Vec3 operator*(const float f, const Vec3 v)
{
    return Vec3{ f * v.x, f * v.y, f * v.z };
}


enum class ObjectType { simple, miroir, glass };

struct Ray
{
    Vec3 origin;
    Vec3 direction;

    Vec3 getXIntersect(float& t) const
    {
        return origin + t * direction; //X = O + t * D
    }
};

struct Sphere
{
    Vec3 centre;
    float radius;
    float albedo;
    ObjectType type;

    float intersect(const Ray& r) const
    {
        Vec3 op = centre - r.origin;
        float t;
        double eps = 1e-4;
        float b = op.scalarProduct(r.direction);
        float det = b * b - op.scalarProduct(op) + radius * radius;
        if (det < 0)
            return 0;
        else
            det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};

struct LightSource
{
    Vec3 origin;
    Vec3d color;
    int emission;
};

struct Box
{
    Vec3 TopLeft;
    Vec3 BottomRight;
    vector<Box> children;
    vector<Sphere> objects;
};


double distanceSquared(Vec3 p1, Vec3 p2) {
    return pow((p1.x + p2.x), 2) + pow((p1.y + p2.y), 2) + pow((p1.z + p2.z), 2);
}

double clamp255(double x) { return x < 0 ? 0 : x > 255 ? 255 : x; }
double clamp(double x, double max) {
    if (x <= 0) {
        return 0;
    }
    else {
        return (x * 255) / max;
    }
}


double getLightColor(float t, Ray ray, LightSource lightsource, Sphere s) {
    const double M_PI = 3.14;

    Vec3 X = ray.getXIntersect(t);
    Vec3 d = lightsource.origin - X;
    float d2 = d.normSquared();
    Vec3 w0 = d / sqrt(d2);

    Vec3 n = X - s.centre;
    float n2 = n.normSquared();
    Vec3 N = n / sqrt(n2);

    double f = w0.scalarProduct(N) / M_PI;

    return (lightsource.emission * f * s.albedo) / d2;
}

vector<Sphere> creationSpheres(float width, float height) {
    //Components Object
    float backgroundDepth = 1600;
    vector<Sphere> spheres;
    spheres.push_back(Sphere{ Vec3{ 500, 350, 200 }, 100, 1, ObjectType::miroir }); //Sphere { center, radius, albedo }
    spheres.push_back(Sphere{ Vec3{ 500, 650, 200 }, 100, 1, ObjectType::simple }); //Sphere { center, radius, albedo }
    spheres.push_back(Sphere{ Vec3{ -height - 200, width / 2, backgroundDepth - 600 }, 1000, 1, ObjectType::simple }); //Sphere { center, radius, albedo }
    spheres.push_back(Sphere{ Vec3{ height + height + 200, width / 2, backgroundDepth - 600 }, 1000, 1, ObjectType::simple }); //Sphere { center, radius, albedo }
    spheres.push_back(Sphere{ Vec3{ height / 2, -width + width / 10, backgroundDepth - 600 }, 1000, 1, ObjectType::simple }); //Sphere { center, radius, albedo }
    spheres.push_back(Sphere{ Vec3{ height / 2, width + width - width / 10, backgroundDepth - 600 }, 1000, 1, ObjectType::simple }); //Sphere { center, radius, albedo }
    spheres.push_back(Sphere{ Vec3{ height / 2, width / 2, backgroundDepth - 100 }, 1000, 1, ObjectType::simple }); //Sphere { center, radius, albedo }
    return spheres;
}


vector<Sphere> creationSpheresBox(float width, float height) {
    //Components Object
    float depth = 400;
    vector<Sphere> spheres;
    int nbLines = 10;
    int nbCols = 10;
    float radius = 25;
    for (int i = 0; i < nbLines; i++) {
        for (int j = 0; j < nbCols; j++) {
            spheres.push_back(Sphere{ Vec3{ 100 + (radius * 2.5f) * i, 100 + (radius * 2.5f) * j, depth }, radius , 1 }); //Sphere { center, radius, albedo }
        }
    }
    return spheres;
}


double rayTracingGetLight(const Ray& ray, const Ray& ray2, float d2, float& t_min, const Sphere& sphere, const vector<Sphere>& spheres, const LightSource& lightsource, double& max_Lc) {
    bool visible = true;
    for (const Sphere& s : spheres) {
        float t = s.intersect(ray2);
        if (t > 0 && t < d2) { // Intersection found
            visible = false;
            break;
        }
    }
    if (visible) {
        double Lc = getLightColor(t_min, ray, lightsource, sphere);
        if (max_Lc < Lc)
        {
            max_Lc = Lc;
        }
        return Lc;
    }
    else {
        return 0;
    }
}


double rayTracingRecursive(const Ray& ray, float& offset, const vector<Sphere>& spheres, const LightSource& lightsource, double& max_Lc, int step) {
    if (step > 10) {
        return 0;
    }

    //Intersection Rayon - Sphere
    float t_min = 0;
    Sphere sphere;
    for (const Sphere& s : spheres) {
        float t = s.intersect(ray);
        if (t > 0) { // Intersection found
            if (t_min == 0 || t < t_min) {
                t_min = t;
                sphere = s;
            }
        }
    }

    //Intersection Rayon - Light
    if (t_min > 0) {
        Vec3 X = ray.getXIntersect(t_min);

        if (sphere.type == ObjectType::simple) {
            Vec3 d = lightsource.origin - X;
            float d2 = d.normSquared();
            Vec3 w1 = d / sqrt(d2);
            Ray ray2{ Vec3{ X + offset * w1 }, w1 }; //Ray { origin, direction }
            return rayTracingGetLight(ray, ray2, d2, t_min, sphere, spheres, lightsource, max_Lc);
        }
        else if (sphere.type == ObjectType::miroir) {
            Vec3 N = X.getNormal(sphere.centre);
            Vec3 R = ray.direction - 2 * ray.direction.scalarProduct(N) * N;
            Ray ray2{ Vec3{ X + offset * R }, R }; //Ray { origin, direction }
            step = step + 1;
            return rayTracingRecursive(ray2, offset, spheres, lightsource, max_Lc, step);
        }
        else if (sphere.type == ObjectType::glass) {
            Vec3 N = X.getNormal(sphere.centre);
            Vec3 R = ray.direction - 2 * ray.direction.scalarProduct(N) * N;
            Ray ray2{ Vec3{ X + offset * R }, R }; //Ray { origin, direction }
            step = step + 1;
            double Lc1 = rayTracingRecursive(ray2, offset, spheres, lightsource, max_Lc, step); // reflection


            float eta = 1.33f;
            float R0;
            if (step > 0) {
                R0 = pow(((eta - 1.0f) / (eta + 1.0f)), 2);
            }
            else {
                R0 = pow(((1.0f - eta) / (1.0f + eta)), 2);
            }
            float k = R0 + (1 - R0) * pow((1 + ray.direction.scalarProduct(N)), 5);

            //float k = 1 - eta * eta * (1 - N.scalarProduct(ray.direction) * N.scalarProduct(ray.direction));
            if (k < 0) {
                return Lc1 * k;
            }
            else {
            
                R = eta * ray.direction - (eta * N.scalarProduct(ray.direction) + sqrt(k)) * N;

                Ray ray3{ Vec3{ X + offset * Vec3{ 1 - R.x, 1 - R.y, 1 - R.z} }, Vec3{ 1 - R.x, 1 - R.y, 1 - R.z} }; //Ray { origin, direction }
                
                double Lc2 = rayTracingRecursive(ray3, offset, spheres, lightsource, max_Lc, step); // refraction

                //cout << Lc1 * k + Lc2 * (1 - double(k)) << endl;

                return Lc1 * k + Lc2 * (1 - double(k));
            }
        }
    }
    else {
        return -1;
    }
}

vector<double> rayTracing(Mat& image, float width, float height, const vector<Sphere>& spheres, const LightSource& lightsource, double& max_Lc) {
    Vec3 rayOrigin = Vec3{ height / 2, width / 2, -1000 };
    float offset = 0.02;

    vector<double> Lcs;

    for (int x = 0; x < image.rows; x++)
    {
        for (int y = 0; y < image.cols; y++)
        {
            //Ray 
            Vec3 pixelOrigin = Vec3{ (float)x,(float)y, 0 };
            Vec3 rayDirectionNormalise = pixelOrigin.getNormal(rayOrigin);
            Ray ray{ rayOrigin, rayDirectionNormalise };

            double Lc = rayTracingRecursive(ray, offset, spheres, lightsource, max_Lc, 0);
            double Lc_clamp = Lc + 0.02;
            Lcs.push_back(Lc_clamp);
        }
    }
    return Lcs;
}

void colorImage(Mat& image, const LightSource& lightsource, const vector<double>& Lcs) {
    //Couleur pour l'image
    Vec3b defaultColor = { 255, 255, 40 };
    Vec3b color;
    for (int x = 0; x < image.rows; x++)
    {
        for (int y = 0; y < image.cols; y++)
        {
            int i = x * image.cols + y;
            double Lc = Lcs.at(i);
            if (Lc < 0) {
                //std::cout << Lc << "\n";
                color = defaultColor;
            }
            else {
                for (int i = 0; i < 3; i++) {
                    double intensity = clamp255(Lc * lightsource.color[i]);
                    color[i] = unsigned char(intensity);
                }
            }
            image.at<Vec3b>(x, y) = color;
        }
    }
}


bool calculBox(const vector<Sphere>& spheres, Vec3& minMinimum, Vec3& maxMaximum) {
    minMinimum = { INFINITY, INFINITY, INFINITY };
    maxMaximum = { -INFINITY, -INFINITY, -INFINITY };
    if (spheres.size() > 0) {
        for (const Sphere& s : spheres) {
            Vec3 minimum = Vec3{ s.centre.x - s.radius, s.centre.y - s.radius, s.centre.z - s.radius };
            Vec3 maximum = Vec3{ s.centre.x + s.radius, s.centre.y + s.radius, s.centre.z + s.radius };
            if (minimum.x < minMinimum.x) {
                minMinimum.x = minimum.x;
            }
            if (minimum.y < minMinimum.y) {
                minMinimum.y = minimum.y;
            }
            if (minimum.z < minMinimum.z) {
                minMinimum.z = minimum.z;
            }

            if (maximum.x > maxMaximum.x) {
                maxMaximum.x = maximum.x;
            }
            if (maximum.y > maxMaximum.y) {
                maxMaximum.y = maximum.y;
            }
            if (maximum.z > maxMaximum.z) {
                maxMaximum.z = maximum.z;
            }
        }
        return true;
    }
    else {
        return false;
    }
}

void calculBoxObject(const vector<Sphere>& spheres, Box& box) {
    Vec3 minMinimum;
    Vec3 maxMaximum;
    bool res = calculBox(spheres, minMinimum, maxMaximum);
    if (res) {
        box.TopLeft = minMinimum;
        box.BottomRight = maxMaximum;

        Vec3 minimum1 = box.TopLeft;
        Vec3 maximum1 = Vec3{ box.BottomRight.x, box.TopLeft.y + (abs(box.BottomRight.y - box.TopLeft.y)) / 2, box.BottomRight.z };
        Vec3 minimum2 = Vec3{ box.TopLeft.x, box.TopLeft.y + (abs(box.BottomRight.y - box.TopLeft.y)) / 2, box.TopLeft.z };
        Vec3 maximum2 = box.BottomRight;

        Box b1 = Box{ minimum1, maximum1, vector<Box>{} };
        Box b2 = Box{ minimum2, maximum2, vector<Box>{} };

        for (const Sphere& s : spheres) {
            Vec3 minimum = Vec3{ s.centre.x - s.radius, s.centre.y - s.radius, s.centre.z - s.radius };
            Vec3 maximum = Vec3{ s.centre.x + s.radius, s.centre.y + s.radius, s.centre.z + s.radius };
            if (b1.TopLeft.x <= minimum.x && b1.TopLeft.y <= minimum.y && b1.TopLeft.z <= minimum.z && b1.BottomRight.x >= maximum.x && b1.BottomRight.y >= maximum.y && b1.BottomRight.z >= maximum.z) {
                b1.objects.push_back(s);
            }
            if (b2.TopLeft.x <= minimum.x && b2.TopLeft.y <= minimum.y && b2.TopLeft.z <= minimum.z && b2.BottomRight.x >= maximum.x && b2.BottomRight.y >= maximum.y && b2.BottomRight.z >= maximum.z) {
                b2.objects.push_back(s);
            }
        }

        box.children.push_back(b1);
        box.children.push_back(b2);

    }
}


void drawBox(Mat& image, const Box& b) {
    float width = abs(b.TopLeft.y - b.BottomRight.y);
    float height = abs(b.TopLeft.x - b.BottomRight.x);
    //cout << b.TopLeft.y << " ! " << b.TopLeft.x << " ! " << width << " ! " << height << endl;
    cv::rectangle(image, Rect(b.TopLeft.y, b.TopLeft.x, width, height), Scalar(255, 255, 0));
}

void drawBoxes(Mat& image, const Box& box, int i) {
    drawBox(image, box);
    for (const Box& b : box.children) {
        drawBoxes(image, b, i + 1);
    }
    for (const Sphere& s : box.objects) {
        cv::circle(image, Point2f(s.centre.x, s.centre.x), 5, Scalar(0, 255, 255), 1);
    }
}

int main()
{
    //Components Image
    const int width = 1000;
    const int height = 800;
    Mat image = Mat::zeros(height, width, CV_8UC3);

    vector<Sphere> spheres = creationSpheres(width, height);
    LightSource lightsource = { Vec3{ 250, 500, 100}, {100, 150, 150}, 1000000 }; //LightSource { origin, couleur, emission }

    double max_Lc = 0;
    vector<double> Lcs = rayTracing(image, width, height, spheres, lightsource, max_Lc);
    colorImage(image, lightsource, Lcs);

    /*
    Box box;
    calculBoxObject(spheres, box);
    drawBoxes(image, box, 1);
    */

    cv::circle(image, Point2f(lightsource.origin.y, lightsource.origin.x), 3, Scalar(0, 255, 255), 1);

    cv::imwrite("./Images/MyImage.png", image);
    cv::imshow("Display Window", image);

    cv::waitKey(0);
    return 0;
}
