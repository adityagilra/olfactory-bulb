% File: readfile_pairstats_spiketime.m
% Function: Pair and non-pair odour response spectrum comparison - read t
% files version
% Files required: List of units & list of pairs/non-pairs
% Author: Ashesh Dhawale
% Date: 04/10/2009

paircorr = [];
npaircorr = [];

spaircorr = [];
snpaircorr = [];

apaircorr = [];

bin = 5;   % ms
wind = 50; % ms, half window duration


clear odourvect1; clear odourvect2; clear ors1; clear ors2;

% Read in unit list file
[sernum unit odname odnum rep phnum odstart odend trend path] = textread('unitfile_100610.txt', '%u %u %s %u %u %u %u %u %u %q', 'headerlines', 1);

% Read in pair/nonpair list
[u1 u2 rel] = textread('pairlist_100610.txt', '%u %u %s', 'headerlines', 1);

fclose all;

pcount = 1; ncount = 1;

currdir = pwd;


for x = 1 : length(u1)
    
    cd (path{u1(x)});
    t1 = loadtfile (odname{u1(x)}, unit(u1(x)));
    spikeset1 = sortrials (t1 , 3, odnum(u1(x)), rep(u1(x)));
    th1 = psth (spikeset1, 25, trend(u1(x)) * 10000, 0);
    phaset1 = sortphase (spikeset1, odname{u1(x)}, phnum(u1(x)), 3, odnum(u1(x)), rep(u1(x)));
    odourvect1 = resp (spikeset1, phaset1, phnum(u1(x)), 3, odstart(u1(x)), odend(u1(x)));
   
    cd (path{u2(x)});
    t2 = loadtfile (odname{u2(x)}, unit(u2(x)));
    spikeset2 = sortrials (t2 , 3, odnum(u2(x)), rep(u2(x)));
    th2 = psth (spikeset2, 25, trend(u2(x)) * 10000, 0);
    phaset2 = sortphase (spikeset2, odname{u2(x)}, phnum(u2(x)), 3, odnum(u2(x)), rep(u2(x)));
    odourvect2 = resp (spikeset2, phaset2, phnum(u2(x)), 3, odstart(u2(x)), odend(u2(x)));    
    
    odtime = odend(u1(x)) - odstart(u1(x));
    airtime = odstart(u1(x));
    
    stc = zeros(8, ((wind/bin)*2) + 1);
    
    bound = bin * (length(stc(1, :)) - 1) * 10/2;
    
    for od = 1 : odnum(u1(x))
        for r = 1 : rep(u1(x))
            
            spikes1 = spikeset1 {od, r};
            spikes2 = spikeset2 {od, r};
            
            % Unit 1 Auto Air
            
            sp1 = spikes1(find((spikes1 > bound) & (spikes1 < (odstart(u1(x))*10000 - bound))));
            
            for s = 1 : length(sp1)
                
                spx = spikes1(find((spikes1 > (sp1(s) - bound)) & (spikes1 < (sp1(s) + bound)))) - sp1(s);
                shist = hist(spx(spx ~= 0), [-bound : 10 * bin : bound]);
                
                stc(1, :) = stc(1, :) + shist;
                
            end
            
            
            % Unit 2 Auto Air
            
            sp2 = spikes2(find((spikes2 > bound) & (spikes2 < (odstart(u2(x))*10000 - bound))));
            
            for s = 1 : length(sp2)
                
                spx = spikes2(find((spikes2 > (sp2(s) - bound)) & (spikes2 < (sp2(s) + bound)))) - sp2(s);
                shist = hist(spx(spx ~= 0), [-bound : 10 * bin : bound]);
                
                stc(2, :) = stc(2, :) + shist;
                
            end
            
            % Unit 1 Auto Odour
            
            sp1 = spikes1(find((spikes1 > (odstart(u1(x))*10000) + bound) & (spikes1 < (odend(u1(x))*10000 - bound))));
            
            for s = 1 : length(sp1)
                
                spx = spikes1(find((spikes1 > (sp1(s) - bound)) & (spikes1 < (sp1(s) + bound)))) - sp1(s);
                shist = hist(spx(spx ~= 0), [-bound : 10 * bin : bound]);
                
                stc(3, :) = stc(3, :) + shist;
                
            end
                        
            % Unit 2 Auto Odour
            
            sp2 = spikes2(find((spikes2 > (odstart(u2(x))*10000) + bound) & (spikes2 < (odend(u2(x))*10000 - bound))));
            
            for s = 1 : length(sp2)
                
                spx = spikes2(find((spikes2 > (sp2(s) - bound)) & (spikes2 < (sp2(s) + bound)))) - sp2(s);
                shist = hist(spx(spx ~= 0), [-bound : 10 * bin : bound]);
                
                stc(4, :) = stc(4, :) + shist;
                
            end
            
            % Unit 1-2 cross Air
            
            sp1 = spikes1(find((spikes1 > bound) & (spikes1 < (odstart(u1(x))*10000 - bound))));
            sp2 = spikes2(find((spikes2 > bound) & (spikes2 < (odstart(u2(x))*10000 - bound))));
            
            for s = 1 : length(sp1)
                
                spx = spikes2(find((spikes2 > (sp1(s) - bound)) & (spikes2 < (sp1(s) + bound)))) - sp1(s);
                shist = hist(spx, [-bound : 10 * bin : bound]);
                
                stc(5, :) = stc(5, :) + shist;
                
            end
            
            % Unit 2-1 cross Air
            
            sp1 = spikes1(find((spikes1 > bound) & (spikes1 < (odstart(u1(x))*10000 - bound))));
            sp2 = spikes2(find((spikes2 > bound) & (spikes2 < (odstart(u2(x))*10000 - bound))));
            
            for s = 1 : length(sp2)
                
                spx = spikes1(find((spikes1 > (sp2(s) - bound)) & (spikes1 < (sp2(s) + bound)))) - sp2(s);
                shist = hist(spx, [-bound : 10 * bin : bound]);
                
                stc(6, :) = stc(6, :) + shist;
                
            end
            
            % Unit 1-2 cross odour
            
            sp1 = spikes1(find((spikes1 > (odstart(u1(x))*10000) + bound) & (spikes1 < (odend(u1(x))*10000 - bound))));
            sp2 = spikes2(find((spikes2 > (odstart(u2(x))*10000) + bound) & (spikes2 < (odend(u2(x))*10000 - bound))));
            
            for s = 1 : length(sp1)
                
                spx = spikes2(find((spikes2 > (sp1(s) - bound)) & (spikes2 < (sp1(s) + bound)))) - sp1(s);
                shist = hist(spx, [-bound : 10 * bin : bound]);
                
                stc(7, :) = stc(7, :) + shist;
                
            end
            
            % Unit 2-1 cross odour
            
            sp1 = spikes1(find((spikes1 > (odstart(u1(x))*10000) + bound) & (spikes1 < (odend(u1(x))*10000 - bound))));
            sp2 = spikes2(find((spikes2 > (odstart(u2(x))*10000) + bound) & (spikes2 < (odend(u2(x))*10000 - bound))));
            
            for s = 1 : length(sp2)
                
                spx = spikes1(find((spikes1 > (sp2(s) - bound)) & (spikes1 < (sp2(s) + bound)))) - sp2(s);
                shist = hist(spx, [-bound : 10 * bin : bound]);
                
                stc(8, :) = stc(8, :) + shist;
                
            end
            
        end
    end
%     
%     stc(1, (wind/bin) + 1) = 0;
%     stc(2, (wind/bin) + 1) = 0;
%     stc(3, (wind/bin) + 1) = 0;
%     stc(4, (wind/bin) + 1) = 0;
    
%     stc(5, (wind/bin) + 1) = NaN;
%     stc(6, (wind/bin) + 1) = NaN;
%     stc(7, (wind/bin) + 1) = NaN;
%     stc(8, (wind/bin) + 1) = NaN;
    
    switch rel {x}
        
        case 'pair'
            
            paircorr(pcount, :, :) = stc;
            
            pcount = pcount + 1;
            
        case 'npair'
            
            npaircorr(ncount, :, :) = stc;

            ncount = ncount + 1;
    end

end
    
cd(currdir);


sm = 1e-5;  % smoothing param for smoothing spline
sbin = 1;   % bin in ms

for x = 1 : size(paircorr, 1)
    
%     figure(x); 
    
    subplot(2, 2, 1); 
    plot([-wind:bin:wind], squeeze(paircorr(x, 1, :))/sum(squeeze(paircorr(x, 1, :))), 'b');
    hold on;
    plot([-wind:bin:wind], squeeze(paircorr(x, 3, :))/sum(squeeze(paircorr(x, 3, :))), 'r');
    
    pp = csaps([-wind:bin:wind], squeeze(paircorr(x, 1, :))/sum(squeeze(paircorr(x, 1, :))), sm);
    spaircorr(x, 1, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'b');
    pp = csaps([-wind:bin:wind], squeeze(paircorr(x, 3, :))/sum(squeeze(paircorr(x, 3, :))), sm);
    spaircorr(x, 3, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'r');    
    
    
    subplot(2, 2, 2); 
    plot([-wind:bin:wind], squeeze(paircorr(x, 2, :))/sum(squeeze(paircorr(x, 2, :))), 'b');
    hold on;
    plot([-wind:bin:wind], squeeze(paircorr(x, 4, :))/sum(squeeze(paircorr(x, 4, :))), 'r');
    
    pp = csaps([-wind:bin:wind], squeeze(paircorr(x, 2, :))/sum(squeeze(paircorr(x, 2, :))), sm);
    spaircorr(x, 2, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'b');
    pp = csaps([-wind:bin:wind], squeeze(paircorr(x, 4, :))/sum(squeeze(paircorr(x, 4, :))), sm);
    spaircorr(x, 4, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'r');    

    
    subplot(2, 2, 3); 
    plot([-wind:bin:wind], squeeze(paircorr(x, 5, :))/sum(squeeze(paircorr(x, 5, :))), 'b');
    hold on;
    plot([-wind:bin:wind], squeeze(paircorr(x, 7, :))/sum(squeeze(paircorr(x, 7, :))), 'r');
    
    pp = csaps([-wind:bin:wind], squeeze(paircorr(x, 5, :))/sum(squeeze(paircorr(x, 5, :))), sm);
    spaircorr(x, 5, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'b');
    pp = csaps([-wind:bin:wind], squeeze(paircorr(x, 7, :))/sum(squeeze(paircorr(x, 7, :))), sm);
    spaircorr(x, 7, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'r');    

    
    subplot(2, 2, 4); 
    plot([-wind:bin:wind], squeeze(paircorr(x, 6, :))/sum(squeeze(paircorr(x, 6, :))), 'b');
    hold on;
    plot([-wind:bin:wind], squeeze(paircorr(x, 8, :))/sum(squeeze(paircorr(x, 8, :))), 'r');
    
    pp = csaps([-wind:bin:wind], squeeze(paircorr(x, 6, :))/sum(squeeze(paircorr(x, 6, :))), sm);
    spaircorr(x, 6, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'b');
    pp = csaps([-wind:bin:wind], squeeze(paircorr(x, 8, :))/sum(squeeze(paircorr(x, 8, :))), sm);
    spaircorr(x, 8, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'r');    

end


for x = 1 : size(npaircorr, 1)
    
%     figure(size(paircorr, 1) + x); 
    
    subplot(2, 2, 1); 
    plot([-wind:bin:wind], squeeze(npaircorr(x, 1, :))/sum(squeeze(npaircorr(x, 1, :))), 'b');
    hold on;
    plot([-wind:bin:wind], squeeze(npaircorr(x, 3, :))/sum(squeeze(npaircorr(x, 3, :))), 'r');
    
    pp = csaps([-wind:bin:wind], squeeze(npaircorr(x, 1, :))/sum(squeeze(npaircorr(x, 1, :))), sm);
    snpaircorr(x, 1, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'b');
    pp = csaps([-wind:bin:wind], squeeze(npaircorr(x, 3, :))/sum(squeeze(npaircorr(x, 3, :))), sm);
    snpaircorr(x, 3, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'r');    

    
    subplot(2, 2, 2); 
    plot([-wind:bin:wind], squeeze(npaircorr(x, 2, :))/sum(squeeze(npaircorr(x, 2, :))), 'b');
    hold on;
    plot([-wind:bin:wind], squeeze(npaircorr(x, 4, :))/sum(squeeze(npaircorr(x, 4, :))), 'r');
    
    pp = csaps([-wind:bin:wind], squeeze(npaircorr(x, 2, :))/sum(squeeze(npaircorr(x, 2, :))), sm);
    snpaircorr(x, 2, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'b');
    pp = csaps([-wind:bin:wind], squeeze(npaircorr(x, 4, :))/sum(squeeze(npaircorr(x, 4, :))), sm);
    snpaircorr(x, 4, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'r');    

    
    subplot(2, 2, 3); 
    plot([-wind:bin:wind], squeeze(npaircorr(x, 5, :))/sum(squeeze(npaircorr(x, 5, :))), 'b');
    hold on;
    plot([-wind:bin:wind], squeeze(npaircorr(x, 7, :))/sum(squeeze(npaircorr(x, 7, :))), 'r');
    
    pp = csaps([-wind:bin:wind], squeeze(npaircorr(x, 5, :))/sum(squeeze(npaircorr(x, 5, :))), sm);
    snpaircorr(x, 5, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'b');
    pp = csaps([-wind:bin:wind], squeeze(npaircorr(x, 7, :))/sum(squeeze(npaircorr(x, 7, :))), sm);
    snpaircorr(x, 7, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'r');    

    
    subplot(2, 2, 4); 
    plot([-wind:bin:wind], squeeze(npaircorr(x, 6, :))/sum(squeeze(npaircorr(x, 6, :))), 'b');
    hold on;
    plot([-wind:bin:wind], squeeze(npaircorr(x, 8, :))/sum(squeeze(npaircorr(x, 8, :))), 'r');
    
    pp = csaps([-wind:bin:wind], squeeze(npaircorr(x, 6, :))/sum(squeeze(npaircorr(x, 6, :))), sm);
    snpaircorr(x, 6, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'b');
    pp = csaps([-wind:bin:wind], squeeze(npaircorr(x, 8, :))/sum(squeeze(npaircorr(x, 8, :))), sm);
    snpaircorr(x, 8, :) = fnval(pp, [-wind:sbin:wind]);
    plot ([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'r');    

    
end


% Pair Stats

astats = zeros(size(paircorr, 1), 5);
ostats = zeros(size(paircorr, 1), 5);

leave = 50; % edge window in which not to detect peaks

for x = 1: size(spaircorr, 1)
    
    atemp = squeeze(spaircorr(x, 5, :));
    otemp = squeeze(spaircorr(x, 7, :));
    
    % Peak Height and position
    [astats(x, 1) astats(x, 2)] = max(atemp((leave/sbin):(end - (leave/sbin))));
    [ostats(x, 1) ostats(x, 2)] = max(otemp((leave/sbin):(end - (leave/sbin))));
    
    astats(x, 2) = astats(x, 2) + (leave/sbin);
    ostats(x, 2) = ostats(x, 2) + (leave/sbin);
    
    % Valley Height
    astats(x, 3) = mean([min(atemp(1:astats(x, 2))) min(atemp(astats(x, 2):end))]);
    ostats(x, 3) = mean([min(otemp(1:ostats(x, 2))) min(otemp(ostats(x, 2):end))]);
    
    % Peak/Valley Ratio
    
    astats(x, 4) = astats(x, 1) / astats(x, 3);
    ostats(x, 4) = ostats(x, 1) / ostats(x, 3);
    
    % FWHM
    
    noiseseg = squeeze(paircorr(x, 5, :)) / sum(squeeze(paircorr(x, 5, :)));
    noise = noiseseg - decimate(squeeze(spaircorr(x, 5, :)), bin/sbin);
    
    if (astats(x, 1) - astats(x, 3)) > 2 * std(noise)
        hm = ((astats(x, 1) - astats(x, 3)) / 2) + astats(x, 3);
        [crap valead] = min(atemp(1:astats(x, 2)));
        [crap vatail] = min(atemp(astats(x, 2):end));
        leadseg = find(atemp(valead:astats(x, 2)) <= hm);
        tailseg = find(atemp(astats(x, 2):(astats(x, 2) + vatail - 1)) <= hm); 

        astats(x, 5) = ((tailseg(1) + astats(x, 2)) - (leadseg(end) + valead)) * sbin;

        hm = ((ostats(x, 1) - ostats(x, 3)) / 2) + ostats(x, 3);
        [crap valead] = min(otemp(1:ostats(x, 2)));
        [crap vatail] = min(otemp(ostats(x, 2):end));
        leadseg = find(otemp(valead:ostats(x, 2)) <= hm);
        tailseg = find(otemp(ostats(x, 2):(ostats(x, 2) + vatail - 1)) <= hm); 

        ostats(x, 5) = ((tailseg(1) + ostats(x, 2)) - (leadseg(end) + valead)) * sbin;        
    else
        astats(x, :) = [NaN NaN NaN NaN NaN];
        ostats(x, :) = [NaN NaN NaN NaN NaN];
    end
    
end

astats(:, 2) = (astats(:, 2) * sbin) - wind;
ostats(:, 2) = (ostats(:, 2) * sbin) - wind;


% Non-Pair Stats

anstats = zeros(size(npaircorr, 1), 5);
onstats = zeros(size(npaircorr, 1), 5);

leave = 50; % edge window in which not to detect peaks

for x = 1: size(snpaircorr, 1)
    
    atemp = squeeze(snpaircorr(x, 5, :));
    otemp = squeeze(snpaircorr(x, 7, :));
    
    % Peak Height and position
    [anstats(x, 1) anstats(x, 2)] = max(atemp((leave/sbin):(end - (leave/sbin))));
    [onstats(x, 1) onstats(x, 2)] = max(otemp((leave/sbin):(end - (leave/sbin))));
    
    anstats(x, 2) = anstats(x, 2) + (leave/sbin);
    onstats(x, 2) = onstats(x, 2) + (leave/sbin);
    
    % Valley Height
    anstats(x, 3) = mean([min(atemp(1:anstats(x, 2))) min(atemp(anstats(x, 2):end))]);
    onstats(x, 3) = mean([min(otemp(1:onstats(x, 2))) min(otemp(onstats(x, 2):end))]);
    
    % Peak/Valley Ratio
    
    anstats(x, 4) = anstats(x, 1) / anstats(x, 3);
    onstats(x, 4) = onstats(x, 1) / onstats(x, 3);
    
    % FWHM
    
    noiseseg = squeeze(npaircorr(x, 5, :)) / sum(squeeze(npaircorr(x, 5, :)));
    noise = noiseseg - decimate(squeeze(snpaircorr(x, 5, :)), bin/sbin);
    
    if (anstats(x, 1) - anstats(x, 3)) > 2 * std(noise)
        hm = ((anstats(x, 1) - anstats(x, 3)) / 2) + anstats(x, 3);
        [crap valead] = min(atemp(1:anstats(x, 2)));
        [crap vatail] = min(atemp(anstats(x, 2):end));
        leadseg = find(atemp(valead:anstats(x, 2)) <= hm);
        tailseg = find(atemp(anstats(x, 2):(anstats(x, 2) + vatail - 1)) <= hm); 

        anstats(x, 5) = ((tailseg(1) + anstats(x, 2)) - (leadseg(end) + valead)) * sbin;

        hm = ((onstats(x, 1) - onstats(x, 3)) / 2) + onstats(x, 3);
        [crap valead] = min(otemp(1:onstats(x, 2)));
        [crap vatail] = min(otemp(onstats(x, 2):end));
        leadseg = find(otemp(valead:onstats(x, 2)) <= hm);
        tailseg = find(otemp(onstats(x, 2):(onstats(x, 2) + vatail - 1)) <= hm); 

        onstats(x, 5) = ((tailseg(1) + onstats(x, 2)) - (leadseg(end) + valead)) * sbin;        
    else
        anstats(x, :) = [NaN NaN NaN NaN NaN];
        onstats(x, :) = [NaN NaN NaN NaN NaN];
    end
    
end

anstats(:, 2) = (anstats(:, 2) * sbin) - wind;
onstats(:, 2) = (onstats(:, 2) * sbin) - wind;


% Autocorr - get spike time corr for all unique units

for x = 1 : length(sernum)
    
    cd (path{x});
    t1 = loadtfile (odname{x}, unit(x));
    spikeset1 = sortrials (t1 , 3, odnum(x), rep(x));
    th1 = psth (spikeset1, 25, trend(x) * 10000, 0);
    phaset1 = sortphase (spikeset1, odname{x}, phnum(x), 3, odnum(x), rep(x));
    odourvect1 = resp (spikeset1, phaset1, phnum(x), 3, odstart(x), odend(x));
    
    odtime = odend(x) - odstart(x);
    airtime = odstart(x);
    
    stc = zeros(8, ((wind/bin)*2) + 1);
    
    bound = bin * (length(stc(1, :)) - 1) * 10/2;
    
    for od = 1 : odnum(x)
        for r = 1 : rep(x)
            
            spikes1 = spikeset1 {od, r};
                        
            % Unit 1 Auto Air
            
            sp1 = spikes1(find((spikes1 > bound) & (spikes1 < (odstart(x)*10000 - bound))));
            
            for s = 1 : length(sp1)
                
                spx = spikes1(find((spikes1 > (sp1(s) - bound)) & (spikes1 < (sp1(s) + bound)))) - sp1(s);
                shist = hist(spx(spx ~= 0), [-bound : 10 * bin : bound]);
                
                stc(1, :) = stc(1, :) + shist;
                
            end
            
            % Unit 1 Auto Odour
            
            sp1 = spikes1(find((spikes1 > (odstart(x)*10000) + bound) & (spikes1 < (odend(x)*10000 - bound))));
            
            for s = 1 : length(sp1)
                
                spx = spikes1(find((spikes1 > (sp1(s) - bound)) & (spikes1 < (sp1(s) + bound)))) - sp1(s);
                shist = hist(spx(spx ~= 0), [-bound : 10 * bin : bound]);
                
                stc(2, :) = stc(2, :) + shist;
                
            end

        end
    end
 
%     stc(1, (wind/bin) + 1) = 0;
%     stc(2, (wind/bin) + 1) = 0;

    apaircorr(x,:,:) = stc;
    
end

cd(currdir);

sapaircorr = [];
leave = 30; % section of autocorr to not fit in ms

for x = 1 : size(apaircorr, 1)
    
%     figure(size(paircorr, 1) + size(npaircorr, 1) + x); 
    
    plot([-wind:bin:wind], squeeze(apaircorr(x, 1, :))/sum(squeeze(apaircorr(x, 1, :))), 'b');
    hold on;
    plot([-wind:bin:wind], squeeze(apaircorr(x, 2, :))/sum(squeeze(apaircorr(x, 2, :))), 'r');
    
    pp = csaps([[-wind:bin:-leave] [leave:bin:wind]], squeeze(apaircorr(x, 1, [[1:(((wind - leave)/bin) + 1)] [(((wind + leave)/bin) + 1):end]]))/sum(squeeze(apaircorr(x, 1, :))), sm);
    sapaircorr(x, 1, :) = fnval(pp, [-wind:sbin:wind]);
    plot([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'b');
    
    pp = csaps([[-wind:bin:-leave] [leave:bin:wind]], squeeze(apaircorr(x, 2, [[1:(((wind - leave)/bin) + 1)] [(((wind + leave)/bin) + 1):end]]))/sum(squeeze(apaircorr(x, 2, :))), sm);
    sapaircorr(x, 2, :) = fnval(pp, [-wind:sbin:wind]);
    plot([-wind:bin:wind], fnval(pp, [-wind:bin:wind]), 'r');
    hold off;
    
end


% Autocorr Stats

aastats = [];
oastats = [];

leave = 100; % edge window in which not to detect peaks

for x = 1: size(sapaircorr, 1)
    
    atemp = squeeze(sapaircorr(x, 1, :));
    otemp = squeeze(sapaircorr(x, 2, :));
    
    % Peak Height and position
    [aastats(x, 1) aastats(x, 2)] = max(atemp((leave/sbin):(end - (leave/sbin))));
    [oastats(x, 1) oastats(x, 2)] = max(otemp((leave/sbin):(end - (leave/sbin))));
    
    aastats(x, 2) = aastats(x, 2) + (leave/sbin);
    oastats(x, 2) = oastats(x, 2) + (leave/sbin);
    
    % Valley Height
    aastats(x, 3) = mean([min(atemp(1:aastats(x, 2))) min(atemp(aastats(x, 2):end))]);
    oastats(x, 3) = mean([min(otemp(1:oastats(x, 2))) min(otemp(oastats(x, 2):end))]);
    
    % Peak/Valley Ratio
    
    aastats(x, 4) = aastats(x, 1) / aastats(x, 3);
    oastats(x, 4) = oastats(x, 1) / oastats(x, 3);
    
    % FWHM
    
    noiseseg = squeeze(apaircorr(x, 1, :)) / sum(squeeze(apaircorr(x, 1, :)));
    noise = noiseseg - decimate(squeeze(sapaircorr(x, 1, :)), bin/sbin);
    
    if (aastats(x, 1) - aastats(x, 3)) > 2 * std(noise)
        hm = ((aastats(x, 1) - aastats(x, 3)) / 2) + aastats(x, 3);
        [crap valead] = min(atemp(1:aastats(x, 2)));
        [crap vatail] = min(atemp(aastats(x, 2):end));
        leadseg = find(atemp(valead:aastats(x, 2)) <= hm);
        tailseg = find(atemp(aastats(x, 2):(aastats(x, 2) + vatail - 1)) <= hm); 

        aastats(x, 5) = ((tailseg(1) + aastats(x, 2)) - (leadseg(end) + valead)) * sbin;

        hm = ((oastats(x, 1) - oastats(x, 3)) / 2) + oastats(x, 3);
        [crap valead] = min(otemp(1:oastats(x, 2)));
        [crap vatail] = min(otemp(oastats(x, 2):end));
        leadseg = find(otemp(valead:oastats(x, 2)) <= hm);
        tailseg = find(otemp(oastats(x, 2):(oastats(x, 2) + vatail - 1)) <= hm); 
        
        oastats(x, 5) = ((tailseg(1) + oastats(x, 2)) - (leadseg(end) + valead)) * sbin;        
    else
        aastats(x, :) = [NaN NaN NaN NaN NaN];
        oastats(x, :) = [NaN NaN NaN NaN NaN];
    end
    
end

aastats(:, 2) = (aastats(:, 2) * sbin) - wind;
oastats(:, 2) = (oastats(:, 2) * sbin) - wind;

%---------------------------------------------------------------

% Calculate mean curves

meancorr = zeros(6, ((wind/bin)*2) + 1);
stecorr = zeros(6, ((wind/bin)*2) + 1);

aex = [1:40];
pex = [1:20];
nex = [1 3 4 5 6 7 8 9 10 11 12 13 14 15];

% Autocorr
for x = 1 : size(apaircorr, 1)
    apaircorr(x, 1, :) = squeeze(apaircorr(x, 1, :))'/sum(squeeze(apaircorr(x, 1, :)));
    apaircorr(x, 2, :) = squeeze(apaircorr(x, 2, :))'/sum(squeeze(apaircorr(x, 2, :)));
end

meancorr(1, :) = squeeze(mean(apaircorr(aex, 1, :)));
stecorr(1, :) = squeeze(std(apaircorr(aex, 1, :)))/sqrt(length(aex));
meancorr(2, :) = squeeze(mean(apaircorr(aex, 2, :)));
stecorr(2, :) = squeeze(std(apaircorr(aex, 2, :)))/sqrt(length(aex));

% Pair cross-corr
for x = 1 : size(paircorr, 1)
    paircorr(x, 5, :) = squeeze(paircorr(x, 5, :))'/sum(squeeze(paircorr(x, 5, :)));
    paircorr(x, 7, :) = squeeze(paircorr(x, 7, :))'/sum(squeeze(paircorr(x, 7, :)));
end

meancorr(3, :) = squeeze(mean(paircorr(pex, 5, :)));
stecorr(3, :) = squeeze(std(paircorr(pex, 5, :)))/sqrt(length(pex));
meancorr(4, :) = squeeze(mean(paircorr(pex, 7, :)));
stecorr(4, :) = squeeze(std(paircorr(pex, 7, :)))/sqrt(length(pex));

meancorr(3, (wind/bin) + 1) = NaN;
stecorr(3, (wind/bin) + 1) = NaN;
meancorr(4, (wind/bin) + 1) = NaN;
stecorr(4, (wind/bin) + 1) = NaN;

% Non-pair cross-corr
for x = 1 : size(npaircorr, 1)
    npaircorr(x, 5, :) = squeeze(npaircorr(x, 5, :))'/sum(squeeze(npaircorr(x, 5, :)));
    npaircorr(x, 7, :) = squeeze(npaircorr(x, 7, :))'/sum(squeeze(npaircorr(x, 7, :)));
end

meancorr(5, :) = squeeze(mean(npaircorr(nex, 5, :)));
stecorr(5, :) = squeeze(std(npaircorr(nex, 5, :)))/sqrt(length(nex));
meancorr(6, :) = squeeze(mean(npaircorr(nex, 7, :)));
stecorr(6, :) = squeeze(std(npaircorr(nex, 7, :)))/sqrt(length(nex));

meancorr(5, (wind/bin) + 1) = NaN;
stecorr(5, (wind/bin) + 1) = NaN;
meancorr(6, (wind/bin) + 1) = NaN;
stecorr(6, (wind/bin) + 1) = NaN;


% Plot stats figures

close all;

% Peak Stats

figure;
hold on;

% Auto
bar(1, nanmean(aastats(:, 1)), 'k');
errorbar(1, nanmean(aastats(:, 1)), nanstd(aastats(:, 1))/sqrt(sum(~isnan(aastats(:, 1)))), 'k');

bar(2, nanmean(oastats(:, 1)), 'k');
errorbar(2, nanmean(oastats(:, 1)), nanstd(oastats(:, 1))/sqrt(sum(~isnan(oastats(:, 1)))), 'k');

% Sisters
bar(3, nanmean(astats(:, 1)), 'k');
errorbar(3, nanmean(astats(:, 1)), nanstd(astats(:, 1))/sqrt(sum(~isnan(astats(:, 1)))), 'k');

bar(4, nanmean(ostats(:, 1)), 'k');
errorbar(4, nanmean(ostats(:, 1)), nanstd(ostats(:, 1))/sqrt(sum(~isnan(ostats(:, 1)))), 'k');

% Non-sisters
bar(5, nanmean(anstats(:, 1)), 'k');
errorbar(5, nanmean(anstats(:, 1)), nanstd(anstats(:, 1))/sqrt(sum(~isnan(anstats(:, 1)))), 'k');

bar(6, nanmean(onstats(:, 1)), 'k');
errorbar(6, nanmean(onstats(:, 1)), nanstd(onstats(:, 1))/sqrt(sum(~isnan(onstats(:, 1)))), 'k');

hold off;


% Valley Stats

figure;
hold on;

% Auto
bar(1, nanmean(aastats(:, 3)), 'k');
errorbar(1, nanmean(aastats(:, 3)), nanstd(aastats(:, 3))/sqrt(sum(~isnan(aastats(:, 3)))), 'k');

bar(2, nanmean(oastats(:, 3)), 'k');
errorbar(2, nanmean(oastats(:, 3)), nanstd(oastats(:, 3))/sqrt(sum(~isnan(oastats(:, 3)))), 'k');

% Sisters
bar(3, nanmean(astats(:, 3)), 'k');
errorbar(3, nanmean(astats(:, 3)), nanstd(astats(:, 3))/sqrt(sum(~isnan(astats(:, 3)))), 'k');

bar(4, nanmean(ostats(:, 3)), 'k');
errorbar(4, nanmean(ostats(:, 3)), nanstd(ostats(:, 3))/sqrt(sum(~isnan(ostats(:, 3)))), 'k');

% Non-sisters
bar(5, nanmean(anstats(:, 3)), 'k');
errorbar(5, nanmean(anstats(:, 3)), nanstd(anstats(:, 3))/sqrt(sum(~isnan(anstats(:, 3)))), 'k');

bar(6, nanmean(onstats(:, 3)), 'k');
errorbar(6, nanmean(onstats(:, 3)), nanstd(onstats(:, 3))/sqrt(sum(~isnan(onstats(:, 3)))), 'k');

hold off;


% FWHM Stats

figure;
hold on;

% Auto
bar(1, nanmean(aastats(:, 5)), 'k');
errorbar(1, nanmean(aastats(:, 5)), nanstd(aastats(:, 5))/sqrt(sum(~isnan(aastats(:, 5)))), 'k');

bar(2, nanmean(oastats(:, 5)), 'k');
errorbar(2, nanmean(oastats(:, 5)), nanstd(oastats(:, 5))/sqrt(sum(~isnan(oastats(:, 5)))), 'k');

% Sisters
bar(3, nanmean(astats(:, 5)), 'k');
errorbar(3, nanmean(astats(:, 5)), nanstd(astats(:, 5))/sqrt(sum(~isnan(astats(:, 5)))), 'k');

bar(4, nanmean(ostats(:, 5)), 'k');
errorbar(4, nanmean(ostats(:, 5)), nanstd(ostats(:, 5))/sqrt(sum(~isnan(ostats(:, 5)))), 'k');

% Non-sisters
bar(5, nanmean(anstats(:, 5)), 'k');
errorbar(5, nanmean(anstats(:, 5)), nanstd(anstats(:, 5))/sqrt(sum(~isnan(anstats(:, 5)))), 'k');

bar(6, nanmean(onstats(:, 5)), 'k');
errorbar(6, nanmean(onstats(:, 5)), nanstd(onstats(:, 5))/sqrt(sum(~isnan(onstats(:, 5)))), 'k');

hold off;


% Peak Shift Stats

figure;
hold on;

% Sisters
bar(3, nanmean(abs(astats(:, 2))), 'k');
errorbar(3, nanmean(abs(astats(:, 2))), nanstd(abs(astats(:, 2)))/sqrt(sum(~isnan(astats(:, 2)))), 'k');

bar(4, nanmean(abs(ostats(:, 2))), 'k');
errorbar(4, nanmean(abs(ostats(:, 2))), nanstd(abs(ostats(:, 2)))/sqrt(sum(~isnan(ostats(:, 2)))), 'k');

% Non-sisters
bar(5, nanmean(abs(anstats(:, 2))), 'k');
errorbar(5, nanmean(abs(anstats(:, 2))), nanstd(abs(anstats(:, 2)))/sqrt(sum(~isnan(anstats(:, 2)))), 'k');

bar(6, nanmean(abs(onstats(:, 2))), 'k');
errorbar(6, nanmean(abs(onstats(:, 2))), nanstd(abs(onstats(:, 2)))/sqrt(sum(~isnan(onstats(:, 2)))), 'k');