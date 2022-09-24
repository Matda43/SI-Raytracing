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
};

Vec3 operator*(const float f, const Vec3 v)
{
    return Vec3{ f * v.x, f * v.y, f * v.z };
}

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


double distanceSquared(Vec3 p1, Vec3 p2) {
    return pow((p1.x + p2.x), 2) + pow((p1.y + p2.y), 2) + pow((p1.z + p2.z), 2);
}

double clamp255(double x) { return x < 0 ? 0 : x > 255 ? 255 : x; }
double clamp1(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
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

int main()
{
    //Components Image
    Mat image = Mat::zeros(800, 1000, CV_8UC3);
    Vec3b defaultColor = { 40, 40, 40 };
    Vec3b color;

    //Components Object
    vector<Sphere> spheres;
    spheres.push_back(Sphere{ Vec3{ 300, 600, 400 }, 100, 1 }); //Sphere { center, radius, albedo }
    spheres.push_back(Sphere{ Vec3{ 500, 400, 600 }, 150, 1 }); //Sphere { center, radius, albedo }
    LightSource lightsource = { Vec3{ 200, 700, 100}, {100, 150, 150}, 1000000 }; //LightSource { origin, couleur, emission }

    vector<double> Lcs;
    double max_Lc = 0;

    for (int x = 0; x < image.rows; x++)
    {
        for (int y = 0; y < image.cols; y++)
        {
            //Ray 
            Ray ray{ Vec3{ (float)x, (float)y, 0 }, Vec3{ 0, 0, 1 } }; //Ray { origin, direction }

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

            if (t_min > 0) {
                Vec3 X = ray.getXIntersect(t_min);
                Vec3 d = lightsource.origin - X;
                float d2 = d.normSquared();
                Vec3 w1 = d / sqrt(d2);

                Ray ray2{ Vec3{ (float)X.x, (float)X.y, X.z }, w1 }; //Ray { origin, direction }

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
                    Lcs.push_back(Lc);
                }
                else {
                    Lcs.push_back(0);
                }
            }
            else {
                Lcs.push_back(-1);
            }
        }
    }

    cout << max_Lc << endl;

    for (int x = 0; x < image.rows; x++)
    {
        for (int y = 0; y < image.cols; y++)
        {
            int i = x * image.cols + y;
            double Lc = Lcs.at(i);
            if (Lc < 0) {
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

    cv::circle(image, Point2f(lightsource.origin.y, lightsource.origin.x), 3, Scalar(0, 255, 255), 1);


    cv::imwrite("./Images/Shadow.png", image);
    cv::imshow("Display Window", image);

    cv::waitKey(0);
    return 0;
}