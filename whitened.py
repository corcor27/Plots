#! /usr/bin/env python
import h5py
import copy
import argparse
from pycbc import DYN_RANGE_FAC, pnutils
from pycbc.io import InferenceFile
from pycbc.inference import option_utils
from pycbc.types import TimeSeries, FrequencySeries
from pycbc import strain as pystrain
from pycbc import waveform
from pycbc.filter import highpass_fir, matched_filter
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from pycbc.filter import match


def get_ylim(data, times, tmin, tmax):
    selected = data[(times >= tmin) & (times < tmax)]
    return 1.1*selected.min(), 1.1*selected.max()


parser = argparse.ArgumentParser()
option_utils.add_inference_results_option_group(parser)
parser.add_argument('--injection-file')
parser.add_argument('--plot-map-waveforms', action='store_true', default=False)
parser.add_argument('--output-file', required=True)
parser.add_argument('--min-xlim')
opts = parser.parse_args()
fp, parameters, _, samples = option_utils.results_from_cli(opts)

## Extract file location from argument
directorys=opts.output_file
folder=directorys[:21]
print folder

## Load dictionary and make array of injected parameters
dic_name="paramDict.npy"
dict_load=folder+dic_name
injected=numpy.load("%s" % dict_load).item()
## Generate array manually
inj_vals=list([1126259462.0,     # time
              injected["mass1"], # mass1
              injected["mass2"], # mass2
              injected["spin1_a"], # spin 1 magnitude
              0, # spin1 azimuthal
              injected["spin1_polar"], # spin 1 polar
              injected["spin2_a"], # spin 2 magnitude
              0, # spin2 azimuthal
              injected["spin2_polar"], # spin 2 polar
              injected["distance"], # distance
              1.5, # coa_phase]
              injected["inclination"], # inclination
              injected["polarization"], # polarisation
              injected["ra"],
              injected["dec"]])

snr_list=[]

fig = pyplot.figure()
ii = 0
colors = {'H1': 'r', 'L1': 'g'}
for ifo in ['H1', 'L1']:
    ii += 1
    ax = fig.add_subplot(int('21{}'.format(ii)))
    print ifo
    # get the psd
    print "loading psd"
    psd = FrequencySeries(fp['{}/psds/0'.format(ifo)][:],
        delta_f=fp['{}/psds/0'.format(ifo)].attrs['delta_f'])\
        /DYN_RANGE_FAC**2

    asd = FrequencySeries(numpy.sqrt(psd.numpy()), delta_f=psd.delta_f)
    # get the strain
    print "loading strain"
    stilde = FrequencySeries(fp['{}/stilde'.format(ifo)][:],
        delta_f=fp['{}/stilde'.format(ifo)].attrs['delta_f'],
        epoch=fp['{}/stilde'.format(ifo)].attrs['epoch'])

    print "whitening"
    wh_stilde = FrequencySeries(stilde / asd, delta_f=stilde.delta_f,
                                 epoch=stilde.epoch)
    wh_strain = wh_stilde.to_timeseries()

    # get the MAP values
    print "loading MAP values"
    llrs = fp.read_likelihood_stats(iteration=opts.iteration,
                                    thin_start=opts.thin_start,
                                    thin_end=opts.thin_end,
                                    thin_interval=opts.thin_interval)
    map_idx = (llrs.loglr + llrs.prior).argmax()
    map_values = samples[map_idx]
    varargs = fp.variable_args
    sargs = fp.static_args
    mapvals = [map_values[arg] for arg in varargs]
    print "generating map waveforms"
    genclass = waveform.select_waveform_generator(fp.static_args['approximant'])
    gen = waveform.FDomainDetFrameGenerator(
        genclass,
        detectors=['H1', 'L1'], epoch=stilde.epoch,
        variable_args=varargs,
        **sargs)
    fs = gen.generate(*map(float, mapvals))[ifo]
    if len(fs) < len(psd):
        fs.resize(len(psd))
    elif len(psd) < len(fs):
        fs = fs[:len(psd)]
    fs /= asd
    ts = fs.to_timeseries()

    print "generating injected waveforms"

    genclass = waveform.select_waveform_generator(fp.static_args['approximant'])
    gen = waveform.FDomainDetFrameGenerator(
        genclass,
        detectors=['H1', 'L1'], epoch=stilde.epoch,
        variable_args=varargs,
        **sargs)
    fis = gen.generate(*map(float, inj_vals))[ifo]
    if len(fis) < len(psd):
        fis.resize(len(psd))
    elif len(psd) < len(fis):
        fis = fis[:len(psd)]
    fis /= asd
    ins = fis.to_timeseries()

    print "plotting"

    try:
        gps_time = sargs['tc']
    except KeyError:
        gps_time = map_values['tc']
    xmin = float(opts.min_xlim)
    xmax = 0.05

    # whitened strain
    x = wh_strain.sample_times.numpy()-gps_time
    y = wh_strain
    ax.plot(x, y, colors[ifo], lw=1.5, zorder=1)
    ylim = get_ylim(y, x, xmin, xmax)

    if opts.plot_map_waveforms:
        ax.plot(ts.sample_times.numpy()-gps_time, ts.data, 'k', lw=2, zorder=2)

    if opts.injection_file:
        # get the injection values
        from pycbc import inject, io
        injf = inject.InjectionSet(opts.injection_file)
        ti = injf.make_strain_from_inj_object(injf.table[0], wh_strain.delta_t, ifo, f_lower=gen.current_params['f_lower'])
        fi = ti.to_frequencyseries(delta_f=gen.current_params['delta_f'])
        if len(fi) < len(psd):
            fi.resize(len(psd))
        elif len(psd) < len(fi):
            fi = fi[:len(psd)]
        fi /= asd
        ti = fi.to_timeseries()
        ax.plot(ti.sample_times.numpy()-gps_time, ti.data, 'b-', lw=2, zorder=2)
    m, i = match(ti, ts, psd=psd, low_frequency_cutoff=gen.current_params['f_lower'])
    print "Match between map and injected is %.2f" % m
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ylim)
    ax.text(xmin, ylim[0], 'Match=%.2f' % m, ha='left', va='bottom', fontsize=12)
    ax.set_ylabel('{} whitened strain'.format(ifo))
    if ii == 2:
        ax.set_xlabel('GPS time - {} (s)'.format(gps_time))

    ## Find and save SNR
    snr = matched_filter(ti, y, psd=psd, low_frequency_cutoff=20.0) ## Should be SNR between injected waveform and whitened data
    snr_map=matched_filter(ts, y, psd=psd, low_frequency_cutoff=20.0) ## SNR between MAP waveform and whitened data

    # Remove regions corrupted by filter wraparound
    snr = snr[len(snr) / 4: len(snr) * 3 / 4]
    snr_map = snr[len(snr_map) / 4: len(snr_map) * 3 / 4]
    snr_list.append(snr)
    snr_list.append(snr_map)
'''
savename=opts.output_file
print savename
figname=savename[:-2]+"_SNR"
## Plot and save SNRs
print "saving SNR figures"
jj=0
pyplot.figure()
for ifo in ['H1', 'L1']: ### Will need to manually extend this for Virgo
    jj+=1
    map_snr=snr_list[jj-1]
    inj_snr=snr_list[jj]
    pyplot.subplot(2,1,jj)
    pyplot.plot(map_snr.sample_times,abs(map_snr))
    pyplot.ylabel("SNR")
pyplot.xlabel("Time")
pyplot.savefig("%s.png" % figname)
'''


print "saving MAP values to dictionary"
MAPDic={}
for aa in range(len(varargs)):
    MAPDic[varargs[aa]]=mapvals[aa]

print "save dictionary"
savename=opts.output_file
print savename
savename=savename[:-2]+"_dic"
numpy.save("%s" %  savename, MAPDic)

fp.close()
fig.savefig(opts.output_file, dpi=200, bbox_inches='tight')
