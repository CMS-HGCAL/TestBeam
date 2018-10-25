/********************************************************************
 * File: waveform.cc
 * ------------------------
 *
 * Description:
 * Tools to analyse raw waveforms.
 *
 *
 * Version:
 * Author: Florian Pitters
 *
 *******************************************************************/

#include "HGCal/Reco/interface/MCP_waveform_reco.h"

using namespace std;



// --------------------------------------------------------
// Defintions
// --------------------------------------------------------

float getBaseline(int peak, short* pulse, int nbinsExcludedLeftOfPeak, int nbinsExcludedRightOfPeak) {
    float sum = 0;
    float cnt = 0;

    // calculate mean of all samples outside of peak region
    if (peak + nbinsExcludedRightOfPeak < 1024 && nbinsExcludedRightOfPeak != -1) {
        for (int i = peak + nbinsExcludedRightOfPeak; i < 1000; i++) {
            sum += pulse[i];
            cnt += 1;
        }
    }

    if (peak - nbinsExcludedLeftOfPeak > 0 && nbinsExcludedLeftOfPeak != -1) {
        for (int i = 10; i < peak - nbinsExcludedLeftOfPeak; i++) {
            sum += pulse[i];
            cnt += 1;
        }
    }

    // option for manual range
    if (nbinsExcludedLeftOfPeak == -1 && nbinsExcludedRightOfPeak == -1) {
        for (int i = 10; i < 1014; i++) {
            if (i < 40) {
                sum += pulse[i];
                cnt += 1;
            }
        }
    }

    if (cnt == 0) {
        return -1;
    }

    return sum / cnt;
}

int substractBaseline(int nsamples, short* pulse, int baseline) {
    for (int i = 0; i < nsamples; i++) {
        pulse[i] = pulse[i] - baseline;
    }

    return 0;
}


int findAbsolutePeak(int nsamples, short* pulse, std::string polarity) {
    int xpeak;
    float ypeak;

    // return error if there is no pulse
    if (nsamples <= 0 || !pulse) {
        return -1;
    }

    // loop over all samples looking for min value
    xpeak = 0;
    ypeak = pulse[20];

    for (int i = 21; i < nsamples - 21; i++) {
        if (polarity == "neg") {
            if (pulse[i] < ypeak) {
                ypeak = pulse[i];
                xpeak = i;
            }
        }
        else if (polarity == "pos") {
            if (pulse[i] > ypeak) {
                ypeak = pulse[i];
                xpeak = i;
            }
        }
        else {
            if (abs(pulse[i]) > abs(ypeak)) {
                ypeak = pulse[i];
                xpeak = i;
            }
        }
    }

    return xpeak;
}


int findRealPeak(int nsamples, short* pulse, int rank, string option) {
    float xmin, ymin, tmp;

    // return error if there is no pulse
    if (nsamples <= 0 || !pulse) {
        return -1;
    }

    short pulse_flt[nsamples];
    for (int i = 0; i < nsamples; i++) {
        pulse_flt[i] = 0;
        tmp = 0;

        // filter with moving average
        if (option == "moving_average") {
            for (int j = -rank; j < rank + 1; j++) {
                tmp += pulse[i + j];
            }
            pulse_flt[i] = tmp / (rank + 1);
        }
        // filter with savitzky golay
        else if (option == "savitzky_golay") {
            pulse_flt[i] = 0;
        }
        // don't filter
        else {
            pulse_flt[i] = pulse[i];
        }
    }

    // loop over all samples looking for min value
    xmin = 0;
    ymin = pulse_flt[10];
    for (int i = 10; i < nsamples - 10; i++) {
        if (ymin > pulse_flt[i]) {
            ymin = pulse_flt[i];
            xmin = i;
        }
    }

    return xmin;
}


int findFirstMinAboveNoise(int nsamples, short* pulse, short noise) {
    float xmin, ymin;

    if (nsamples <= 0 || !pulse) {
        return -1;
    }

    xmin = 0;
    for (int i = 10; i < nsamples - 10; i++) {
        if (abs(pulse[i]) > noise) {
            xmin = i;
            break;
        }
    }

    return xmin;
}


float getNoise(int peak, short* pulse, int nbinsExcludedLeftOfPeak, int nbinsExcludedRightOfPeak) {
    float sum = 0;
    float cnt = 0;
    float var = 0;

    // calculate mean and varriance of all samples outside of peak region
    if (peak + nbinsExcludedRightOfPeak < 1024 && nbinsExcludedRightOfPeak != -1) {
        for (int i = peak + nbinsExcludedRightOfPeak; i < 1000; i++) {
            sum += pulse[i];
            var += pow(pulse[i], 2);
            cnt += 1;
        }
    }

    if (peak - nbinsExcludedLeftOfPeak > 0 && nbinsExcludedLeftOfPeak != -1) {
        for (int i = 10; i < peak - nbinsExcludedLeftOfPeak; i++) {
            sum += pulse[i];
            var += pow(pulse[i], 2);
            cnt += 1;
        }
    }

    // option for manual range
    if (nbinsExcludedLeftOfPeak == -1 && nbinsExcludedRightOfPeak == -1) {
        for (int i = 10; i < 1014; i++) {
            if (i < 50 || i > 900) {
                sum += pulse[i];
                var += pow(pulse[i], 2);
                cnt += 1;
            }
        }
    }

    // break if no cnts
    if (cnt == 0) {
        return -1;
    }

    // calculate variance
    var = (cnt * var - pow(sum, 2)) / pow(cnt, 2);

    return sqrt(var);
}


float getPulseIntegral(int peak, short* pulse, int roi_left, int roi_right, string option) {
    float integral;

    integral = 0;

    // calculate sum of all samples in full waveform
    if (option == "full") {
        for (int i = 100; i < 300; i++) {
            integral += pulse[i];
        }
    }

    // calculate sum of all samples within peak region
    else {
        for (int i = peak - roi_left; i < peak + roi_right; i++) {
            integral += pulse[i];
        }
    }

    return integral;
}



int qualityPulse(int nsamples, short* pulse) {
    int fQuality = 1;

    int ipeak;
    float vpeak, noise;

    // return error if there is no pulse
    if (nsamples < 1 || !pulse) {
        fQuality = 0;
    }

    ipeak = findAbsolutePeak(nsamples, pulse, "pos");
    vpeak = pulse[ipeak];
    noise = getNoise(ipeak, pulse, 30, 70);

    if (vpeak < 1. * noise || vpeak < 50) {
        fQuality = 0;
    }

    return fQuality;
}


//modified by T.Q.
peakValues* analysePeak(int nsamples, short* pulse) {
    int id, peak;
    float amp, charge, base, noise, mu_peakfit, amp_peakfit, tlinear15, tlinear30, tlinear45, tlinear60;

    peakValues* ret = new peakValues();

    int roi_left = 1;
    int roi_right = 1;
    float x[roi_left + roi_right + 1], y[roi_left + roi_right + 1], errX[roi_left + roi_right + 1], errY[roi_left + roi_right + 1];

    peak = findAbsolutePeak(nsamples, pulse, "pos");
    amp = pulse[peak];
    base = getBaseline(peak, pulse, 90, -1);
    noise = getNoise(peak, pulse, 90, -1);
    charge = getPulseIntegral(peak, pulse, roi_left, roi_right, "");

    if (qualityPulse(nsamples, pulse) > 0) {
        //gaussian fitting
        for (int i = 0; i <= roi_left + roi_right; i++) {
            x[i] = i - roi_left; //factor 0.2 removed
            y[i] = pulse[peak + i - roi_left];
            errX[i] = 0.0;
            errY[i] = 1.;
        }
        TGraphErrors* gr_pulse = new TGraphErrors(roi_left + roi_right + 1, x, y, errX, errY);

        TF1* fpeak = new TF1("fpeak", "gaus", -roi_left, roi_right);
        fpeak->SetParameter(0, amp);
        gr_pulse->Fit("fpeak", "Q", "", roi_left, roi_right);
        amp_peakfit = fpeak->GetParameter(0);
        mu_peakfit = fpeak->GetParameter(1);
        delete fpeak;

        TF1* flinear = new TF1("flinear", "[0]*x+[1]", -roi_left, 0);
        gr_pulse->Fit("flinear", "Q", "", -roi_left, -1);
        float m_linfit = flinear->GetParameter(0);
        float b_linfit = flinear->GetParameter(1);

        tlinear15 = (0.15*amp_peakfit-b_linfit)/m_linfit;
        tlinear30 = (0.30*amp_peakfit-b_linfit)/m_linfit;
        tlinear45 = (0.45*amp_peakfit-b_linfit)/m_linfit;
        tlinear60 = (0.60*amp_peakfit-b_linfit)/m_linfit;
        ret->fQuality = 1;
    } else {
        amp_peakfit = -1;
        mu_peakfit = 0;
        tlinear15 = -1;
        tlinear30 = -1;
        tlinear45 = -1;
        tlinear60 = -1;
        ret->fQuality = 0;
    }

    ret->peak = peak;      // index of pulse sample holding the largest entry [-]
    ret->amp = amp;        // value of pulse sample holding the largest entry [adc]
    ret->charge = charge;  // integral of either full pulse or roi [adc]
    ret->base = base;      // mean of voltage samples outside of roi [adc]
    ret->noise = noise;    // std of voltage samples outside of roi [adc]
    ret->amppeak = amp_peakfit;    // rise time of peak
    ret->tpeak = mu_peakfit+peak;    // mean of gauss fit around peak
    ret->tlinear15 = tlinear15;
    ret->tlinear30 = tlinear30;
    ret->tlinear45 = tlinear45;
    ret->tlinear60 = tlinear60;
    // ret->quality = fQuality;


    return ret;
}
