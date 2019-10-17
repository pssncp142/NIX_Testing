from NIX_Testing import NIX_Image_List, NIX_Image
from pylab import *
from astropy.io import fits
from lmfit.models import PolynomialModel

data_dir = '/home/ydallilar/Documents/NIX/nixDetBackup'

test_ids = ['PER-124-04-1']
config = {'test_ids' : test_ids, 'data_dir' : data_dir}
NIX_GL = NIX_Image_List(config)

sz = len(NIX_GL)/2-5

ims = zeros([2048, 2048, sz/2])
diffs = zeros([2048, 2048, sz/2])

BP_Mask = NIX_Image('../../BP_mask.fits')

for i in range(sz/2):
    ims[:,:,i] = NIX_GL[2*i+40].getImage(mask=BP_Mask)-NIX_GL[2*i].getImage(mask=BP_Mask)
    diffs[:,:,i] = NIX_GL[2*i+41].getImage(mask=BP_Mask)-NIX_GL[2*i+40].getImage(mask=BP_Mask)

ims[where(ims == 0)] = NaN
diffs[where(diffs == 0)] = NaN

vars = []
signals = []
true_signal = []

exp = (arange(20)+1)*2

for i in range(20):
    for j in range(20):
        for k in range(sz/2):
            true_signal.append(nanmedian(ims[500+50*i-25:500+50*i+25,
                            500+50*j-25:500+50*j+25, 0]*(k+1)))
            vars.append(nanvar(diffs[500+50*i-25:500+50*i+25,
                            500+50*j-25:500+50*j+25, k])/2.)
            signals.append(nanmedian(ims[500+50*i-25:500+50*i+25,
                            500+50*j-25:500+50*j+25, k]))


signals = array(signals)
vars = array(vars)
true_signal = array(true_signal)

order = 5

ndx = where(signals < 16000)[0]

model = PolynomialModel(order)
pars = model.guess(true_signal[ndx], x=signals[ndx])
pars['c0'].set(value=0., vary=False)

out = model.fit(true_signal[ndx], pars, x=signals[ndx])

print out.fit_report()

pol = poly1d([out.best_values['c%d' % i] for i in range(order+1)[::-1]])

xx = linspace(0, 17000, 100)

plot(signals, signals/true_signal, 'k.', markersize=1.)
#plot(signals, true_signal, '.', markersize=1.)
#plot(xx, pol(xx), 'r-')
show()

for i in range(sz/2):
    ims[:,:,i] = pol(NIX_GL[2*i+40].getImage(mask=BP_Mask))-pol(NIX_GL[2*i].getImage(mask=BP_Mask))
    diffs[:,:,i] = pol(NIX_GL[2*i+41].getImage(mask=BP_Mask))-pol(NIX_GL[2*i+40].getImage(mask=BP_Mask))

ims[where(ims == 0)] = NaN
diffs[where(diffs == 0)] = NaN

vars = []
signals = []
true_signal = []

for i in range(20):
    for j in range(20):
        for k in range(sz/2):
            true_signal.append(nanmedian(ims[500+50*i-25:500+50*i+25,
                            500+50*j-25:500+50*j+25, 0]*(k+1)))
            vars.append(nanvar(diffs[500+50*i-25:500+50*i+25,
                            500+50*j-25:500+50*j+25, k]))
            signals.append(nanmedian(ims[500+50*i-25:500+50*i+25,
                            500+50*j-25:500+50*j+25, k]))

order = 1

signals = array(signals)
vars = array(vars)/2.
true_signal = array(true_signal)

ndx = where(signals < 22000)[0]

model = PolynomialModel(order)
pars = model.guess(signals[ndx], x=vars[ndx])
pars['c0'].set(value=0., vary=False)

out = model.fit(signals[ndx], pars, x=vars[ndx])

print out.fit_report()

xx = linspace(0, 1000, 100)
pol = poly1d([out.best_values['c%d' % i] for i in range(order+1)[::-1]])

plt.figure()
plot(vars, signals/vars, '.', markersize=1.)
#plot(xx, pol(xx), 'r-')

order = 1

log_signal = log10(signals)
log_std = log10(sqrt(vars))

model = PolynomialModel(order)
pars = model.guess(log_std[ndx], x=log_signal[ndx])
pars['c1'].set(value=0.5, vary=False)

out = model.fit(log_std[ndx], pars, x=log_signal[ndx])

print out.fit_report()

xx = linspace(3, 4.4, 100)
pol = poly1d([out.best_values['c%d' % i] for i in range(order+1)[::-1]])


plt.figure()
plot(log_signal, log_std, 'k.', markersize=1.)
plot(xx, pol(xx), 'r-')

show()
