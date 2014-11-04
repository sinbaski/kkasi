clear all
close all
names = {'nordea_bank', 'ericsson_b', 'swedbank_a', 'volvo_b'};
prices = NaN(2e+6, 4);

for c = 1 : length(names)
    [A, B] = get_intra_ret_simple(names{c}, '2013-10-11', '2014-08-21', ...
                                            NaN);
    prices(
end
