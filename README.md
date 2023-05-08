This is a short collection of example on how to use boost library to calculate characteristics of a random sample.
The original source I take from this page : http://www.matrixscience.com/msparser/help/group___shapiro_wilk_source_code.html
To install boost on linux you can simply use this command : 
$sudo apt-get install libboost-all-dev
on Windows do your own research.
To build the example on Linux :
1. cd into each directory
2. mkdir -p build
3. cmake ..
4. make

The eigen example needs eigen library too.



    data = np.log(stacked_data['ask_price'][:data_len]) - np.log(stacked_data['bid_price'][:data_len])
    stat, p = shapiro(data)
    variance = np.var(data)
    hist,bins = np.histogram(data,bins=20)
    hist_max = np.argmax(hist)
    dom = bins[hist_max]
    smooth = np.max(np.diff(np.histogram(data, bins='fd')[0]))
    data2 = np.random.normal(0, 1, len(data))
    corr = np.corrcoef(data, data2)[0, 1]
    cov = np.cov(data, data2)[0, 1] 
    acf = sm.tsa.stattools.acf(data)[1:]
    randomess = np.max(acf)
    
    #volume
    data_len = min(len(stacked_data['bid_volume']),len(stacked_data['ask_volume']))
    data = np.log(stacked_data['ask_volume'][:data_len]) - np.log(stacked_data['bid_volume'][:data_len])
    stat_v, p_v = shapiro(data)
    variance_v = np.var(data)
    hist,bins = np.histogram(data,bins=20)
    hist_max_v = np.argmax(hist)
    dom_v = bins[hist_max]
    smooth_v = np.max(np.diff(np.histogram(data, bins='fd')[0]))
    data2 = np.random.normal(0, 1, len(data))
    corr_v = np.corrcoef(data, data2)[0, 1]
    cov_v = np.cov(data, data2)[0, 1] 
    acf = sm.tsa.stattools.acf(data)[1:]
    randomess_v = np.max(acf)
    close = (float(stacked_data['bid_price'][0]) + float(stacked_data['ask_price'][0]))/2
