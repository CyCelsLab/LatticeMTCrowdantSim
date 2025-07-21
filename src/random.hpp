#pragma once
#include<random>
static std::random_device rd;
static std::mt19937 engine(rd());
static std::uniform_real_distribution<double> random_domain(0.0 , 1.0);
#define random random_domain(engine)