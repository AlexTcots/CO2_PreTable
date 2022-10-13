//
// Created by alext on 11/20/2021.
//

#ifndef FHO_PRETABLE_RANDOM_H
#define FHO_PRETABLE_RANDOM_H


class RanPark {
public:
    RanPark(int);
    RanPark(double);
    ~RanPark() {}
    void reset(double, int, int);
    double uniform();
    double gaussian();

private:
    int seed,save;
    double second;
};



#endif //FHO_PRETABLE_RANDOM_H
