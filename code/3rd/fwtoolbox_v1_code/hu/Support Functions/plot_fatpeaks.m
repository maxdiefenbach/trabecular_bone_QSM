figure;clf;
subplot(121); stem(algoParams.species(2).frequency, algoParams.species(2).relAmps,'ro');
title ('fat spectral model'); grid on; xlabel ('Hz'); ylabel ('amplitude'); axis square;

fprintf ('\n');
disp ('echo spacing is/are (msec): ');
deltaTE = diff(algoParams.te)*1000
disp ('water and -CH2- fat phase angles are (deg): ');
if algoParams.sp_mp==1
    ph_angle = (2*pi*algoParams.species(2).frequency(2)*algoParams.te*(180/pi))
else
    ph_angle = (2*pi*algoParams.species(2).frequency(1)*algoParams.te*(180/pi))
end
ph_angle = abs(ph_angle);
disp ('angle difference (deg): ');
diff (ph_angle)