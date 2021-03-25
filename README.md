# DiffusionMoments

[![Build Status](https://travis-ci.com/SosUts/DiffusionMoments.jl.svg?branch=master)](https://travis-ci.com/SosUts/DiffusionMoments.jl)
[![Coverage](https://codecov.io/gh/SosUts/DiffusionMoments.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SosUts/DiffusionMoments.jl)
[![Coverage](https://coveralls.io/repos/github/SosUts/DiffusionMoments.jl/badge.svg?branch=master)](https://coveralls.io/github/SosUts/DiffusionMoments.jl?branch=master)

## Introduction

DiffusionMoments is a package for Mean Maximal Excursion(MME) and Mean Square Displacement(MSD a.k.a regular moment) analysis of particle trajectory.

## Types of moments

### Regular moments(msd)
#### Ensemble averaged msd(EA-MSD)
$$
\overline{{\bm{r}^2(\tau)}} = \frac{1}{N}\sum_{n=1}^N \left [ \bm{x}_n(\tau) - \bm{x}_n(1) \right ]^2
$$

#### Time averaged msd(TA-MSD)
$$
\left< \bm{r}^2(\tau) \right> = \frac{1}{T-\tau}\sum_{t=1}^{T} \left[ \bm{x}(t+\tau) - \bm{x}(t) \right]^2
$$

#### Ensemble averaged TA-MSD

$$
\left< \overline{{\bm{r}^2(\tau)}} \right> = \frac{1}{T-\tau}\sum_{t=1}^{T} \left[ \bm{x}_n(t+\tau) - \bm{x}_n(t) \right]^2
$$

## Usage
```Julia
using DiffusionMoments
using CSV

df = CSV.read(".\\sample.csv")

# ensemble msd(use all possible time lag $\tau$)
eamsd = ensemble_msd(df, :id, :x, :y)

# ensemble time average msd and time average msd
eatamsd, tamsd = ensemble_time_average_msd(
            df, :id, :x, :y, return_tamsd = true
        )

# ensemble time average msd
eatamsd = ensemble_time_average_msd(
            df, :id, :x, :y, return_tamsd = false
        )

# ensemble msd(use only time lag $\tau = 1$)
eaeamsd = ensemble_tamsd(eamsd)
```