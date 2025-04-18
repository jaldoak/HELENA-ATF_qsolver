%%% Plot q-solver profiles

%%% This is a quick script that will make a set of plots comparing the
%%% output of HATF running with F-input vs q-input.
%%% Change paths and file names as needed.

%Change these paths to the the location of "MHD_scripts-master" and
%the qsolver HELENA+ATF folder
script_path = '~/Documents/MHD_codes/MHD_scripts-master';
hatf_path   = '~/Documents/qsolver/HELENA+ATF';
fname0 = 'F2_example'; %F-in file name
fname1 = 'q_example';  %q-in file name


addpath(script_path)
addpath([script_path,'/HATFMA_scripts'])


%% Plot q profiles (of F2_example and q_example)

hatf0 = readHEL_ATF_DATA([hatf_path,'/plot/',fname0,'/',fname0]);
hatf0_out = readHELENA_ATF([hatf_path,'/output/',fname0]);
hatf1 = readHEL_ATF_DATA([hatf_path,'/plot/',fname1,'/',fname1]);
hatf1_out = readHELENA_ATF([hatf_path,'/output/',fname1]);
A=hatf0_out.A;

%Compare q-profiles 
%(note the q-solver q-output is trivially the same as input.)
%(To actually gauge accuracy, F(Psi) output should be plotted)
figure()
plot(hatf0.qprof(:,1),hatf0.qprof(:,2))
hold on
plot(hatf1.qprof(:,1),hatf1.qprof(:,2),'--')
ylabel('q')
xlabel('\Psi')
title('q-match in HATF from F-in vs q-in')
legend('F-in','q-in')


%Plot FDF
figure()
fspline = spline(hatf1.average(:,1),hatf1.average(:,17),hatf0.prof(:,1));
plot(hatf1.prof(:,1),gradient(hatf0.prof(:,3)*A,hatf0.prof(:,1)))
hold on
plot(hatf1.prof(:,1),fspline,'--')
title('qsolver (DF2)-match from qsolver output')
ylabel("(F^2)'")
xlabel('\Psi')
legend('F-in','q-in')

%Integrate D(F^2) output to get plots of F(Psi)
fq2 = zeros(length(hatf1.average(:,1)),1);
for i = 2:length(hatf1.average(:,1))
    fq2(i) = fq2(i-1)+hatf1.average(i,17)*(hatf1.average(i,1)-hatf1.average(i-1,1));
end
c1 = -fq2(end)+hatf1.average(end,16).^2;
Z1 = cumtrapz(hatf1.average(:,1),hatf1.average(:,17));

%Plot F-match from direct F(q) output  
%[NOTE: this plot will not match.  This is here to demonstrate that integrating is necessary.]
figure()
plot(hatf0.axprof(:,4),hatf0.axprof(:,10))
hold on
plot(hatf1.axprof(:,4),hatf1.axprof(:,10),'--')
title('qsolver f-match from direct output')
ylabel('F')
xlabel('\Psi')
legend('Target','q-solver NR50NP33')

%Plot F-match from integrating FDF (normalised)
figure()
plot(hatf0.axprof(:,4),hatf0.axprof(:,10)/max(hatf0.axprof(:,10)))
hold on
plot(hatf1.average(:,1),(Z1+c1).^0.5 / (max(Z1)+c1).^0.5,'--')
title('qsolver f-match from FDF(2) integration (normalised)')
ylabel('F')
xlabel('\Psi')
legend('Target','q-solver NR50NP33')

%Plot F-match from integrating FDF (un-normalised)
% Note: this plot relies on the calculation of A from F(q), hence it is
% a little inaccurate, especially at low resolution.
figure()
plot(hatf0.axprof(:,4),hatf0.axprof(:,10))
hold on
plot(hatf1.average(:,1),(Z1+c1).^0.5,'--')
title('qsolver f-match from FDF(2) integration (un-normalised)')
ylabel('F')
xlabel('\Psi')
legend('Target','q-solver NR50NP33')







