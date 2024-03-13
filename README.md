# Solar activity forecast
The goal of this project is to examine the structure of the time series of solar activity indices and predict the shape and amplitude of the 25th solar cycle. 

# Data sources
To forecast the solar cycle, *monthly mean total sunspot number (SSN)* data is used (source: [WDC-SILSO](https://www.sidc.be/SILSO/datafiles), Royal Observatory of Belgium, Brussels). Data description:
> Monthly mean total sunspot number obtained by taking a simple arithmetic mean of the daily total sunspot number over all days of each calendar month. Monthly means are available only since 1749 because the original observations compiled by Rudolph Wolf were too sparse before that year. A value of -1 indicates missing value. Data before SC 1 (until 1755-02) was dismissed from the dataset.

CSV file contents:
> * Column 1-2: Gregorian calendar date
> * Column 3: Date in fraction of year.
> * Column 4: Monthly mean total sunspot number.
> * Column 5: Monthly mean standard deviation of the input sunspot numbers.

The need to consider solar cycles separately limits the use of data from February 1755 to November 2019 inclusive, which covers 24 cycles of solar activity. The dates of minima and maxima of solar cycles are taken from [NGDC NOAA](https://www.ngdc.noaa.gov/stp/space-weather/solar-data/solar-indices/sunspot-numbers/cycle-data/table_cycle-dates_maximum-minimum.txt) or [Wikipedia](https://en.wikipedia.org/wiki/List_of_solar_cycles).
The 25th cycle has been in the growth phase since December 2019, which gives opportunity to compare the forecast with the *current observational data* (source: [SWPC-NOAA](https://www.swpc.noaa.gov/products/solar-cycle-progression), ISES Solar Cycle Number Progression).

# Methods
The work used a spectral analysis method based on the transformation of a time series into a multidimensional series, decomposition of a multidimensional series into generalized Fourier-Walsh series, one-step extrapolation of each time series and reverse convolution. 
Multidimensional series are obtained by dividing each of the 24 cycles of solar activity into 16 phases. Since the cycle durations in months are different, when dividing cycles into phases, it is necessary to use a weighted average SSN value for the months falling inside or on the border of each phase.
To decompose multidimensional series into a generalized Fourier-Walsh series, the sympy library and the FWHT method are used. To extrapolate each phase based on the solar cycle number, the following conponents must be sequentionally separated, using [normal probability test](https://en.wikipedia.org/wiki/Normal_probability_plot ) with scipy to check the randomness of the residuals of such separation:
> *   Trend: The increasing or decreasing value in the series.
> *   Seasonality: The repeating short-term cycle in the series.
> *   Noise: The random variation in the series.

# Results
The forecast for the 25th solar cycle indicates that the maximum value of 184.9 will be reached by the smoothed over 8.8 months SSN in September 2025. The forecast correlates with observational data up to January 2024.
