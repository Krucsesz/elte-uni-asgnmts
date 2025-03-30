#include <iostream>
#include "Vector2.h"

int main() {
    Vector2<float> v1(3.0f, 4.0f);
    Vector2<float> v2(1.0f, 2.0f);

    Vector2<float> sum = v1 + v2;
    std::cout << "Összeg: " << sum << std::endl;

    Vector2<float> diff = v1 - v2;
    std::cout << "Különbség: " << diff << std::endl;

    Vector2<float> scaled1 = v1 * 2.0f;
    Vector2<float> scaled2 = 2.0f * v1;
    std::cout << "Skalárszorzás (jobbról): " << scaled1 << std::endl;
    std::cout << "Skalárszorzás (balról): " << scaled2 << std::endl;

    Vector2<float> divided = v1 / 2.0f;
    std::cout << "Skalárosztás: " << divided << std::endl;

    float dotProduct = v1.dot(v2);
    std::cout << "Skaláris szorzat (dot product): " << dotProduct << std::endl;

    float sqlength = v1.sqlength();
    std::cout << "v1 hosszának négyzete: " << sqlength << std::endl;

    float length = v1.length();
    std::cout << "v1 hossza: " << length << std::endl;

    Vector2<float> normalized = v1.normalized();
    std::cout << "v1 normalizálva: " << normalized << std::endl;

    v1 += v2;
    std::cout << "v1 += v2: " << v1 << std::endl;

    v1 -= v2;
    std::cout << "v1 -= v2: " << v1 << std::endl;

    v1 *= 2.0f;
    std::cout << "v1 *= 2: " << v1 << std::endl;

    v1 /= 2.0f;
    std::cout << "v1 /= 2: " << v1 << std::endl;

    return 0;
}