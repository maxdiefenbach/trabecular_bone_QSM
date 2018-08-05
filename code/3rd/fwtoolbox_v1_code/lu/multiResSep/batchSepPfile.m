TE = [2604 4356 4356*2-2604]/1000;

fn = '/home_cuneo/recon/data/2008-06-12--Rt-Leg-DOB-10231973/P06144.7';
NS = 128;

W = []; F = []; PSI = [];

for ctr = 1:NS,
    disp(['Process slice ', num2str(ctr)]);
    dt = recon3d(fn, ctr);
    [w,f,p] = multiResSep(TE, dt, 0);

    W = cat(3, W, w);
    F = cat(3, F, f);
    PSI = cat(3, PSI, p);

    save SepResult W F PSI;
end

