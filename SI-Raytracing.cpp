#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>

# define M_PI           3.14159265358979323846

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
};

struct Sphere
{
    Vec3 centre;
    float radius;

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


int main()
{
    //Components Image
    Mat image = Mat::zeros(800, 1000, CV_8UC3);
    Vec3b defaultColor = { 40, 40, 40 };
    Vec3b color;

    //Components Object
    Sphere spheres[2]{
        Sphere{ Vec3{ 400, 300, 0 }, 150 }, //Sphere { center, radius }
        Sphere{ Vec3{ 400, 700, 0 }, 150 } //Sphere { center, radius }
    };
    LightSource lightsource = { Vec3{ 100, 500, 400}, {100, 150, 150}, 1000000 }; //LightSource { origin, emission }

    //Parameters
    int albedo = 1;

    for (int x = 0; x < image.rows; x++)
    {
        for (int y = 0; y < image.cols; y++)
        {
            //Default color
            color = defaultColor;

            //Ray 
            Ray ray{ Vec3{ (float)x, (float)y, 10 }, Vec3{ 0, 0, 1 } }; //Ray { origin, direction }

            //Intersection Rayon - Sphere
            for (Sphere s : spheres) {
                float t = s.intersect(ray);
                if (t > 0) { // Intersection found

                    Vec3 X = ray.origin + t * ray.direction; //X = O + t * D

                    Vec3 d = lightsource.origin - X;
                    float d2 = d.normSquared();
                    Vec3 w0 = d / sqrt(d2);

                    Vec3 n = X - s.centre;
                    float n2 = n.normSquared();
                    Vec3 N = n / sqrt(n2);

                    double f = w0.scalarProduct(N) / M_PI;

                    double Lc = (lightsource.emission * f * albedo) / d2;

                    for (int i = 0; i < 3; i++) {
                        double intensity = clamp255(Lc * lightsource.color[i]);
                        color[i] = unsigned char(intensity);
                    }
                }

            }

            image.at<Vec3b>(x, y) = color;
        }
    }

    cv::circle(image, Point2f(lightsource.origin.y, lightsource.origin.x), 3, Scalar(0, 255, 255), 1);


    cv::imwrite("./Images/2Spheres.png", image);
    cv::imshow("Display Window", image);
    cv::waitKey(0);
    return 0;
}