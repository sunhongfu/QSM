function [par] = readprocpar(fpath)

% READPROCPAR reads in parameter used in FID acquisition
%
%   PAR = READPROCPAR(FPATH) reads the parameters in the the
%   procpar file in the FID directory FPATH
%
% see also WRITEPROCPAR, READFID, WRITEFID


% Written by L Martyn Klassen
% Copyright 2003 Robarts Research Institute

if nargin < 1
   error('READPROCPAR requires one input argument.');
end

% Check for validity of file
fp = fopen([fpath '/procpar'], 'r');
if fp == -1
   fp = fopen([fpath '.fid/procpar'], 'r');
   if fp == -1
      fp = fopen([fpath '.par/procpar'], 'r');
      if fp == -1
         error(['READPROCPAR unable to open ' fpath ' for reading.']);
      else
         fpath = [fpath '.par'];
      end
   else
      fpath = [fpath '.fid'];
   end
end

% Initialize outputs
par.path     = fpath;
par.dp       = [];
par.seqfil   = [];
par.nf       = [];
par.ni       = [];
par.sw       = [];
ar.ns       = [];
par.ne1      = [];
par.ne2      = [];
par.nss       = [];
par.nse       = [];
par.rfspoil  = [];
par.rfphase  = [];
par.nt       = [];
par.np       = [];
par.nv       = [];
par.nv2      = [];
par.nv3      = [];
par.ne       = [];
par.polarity = [];
par.evenecho = [];
par.tr       = [];
par.te       = [];
par.esp      = [];
par.espincr  = [];
par.nprof    = [];
par.tproj    = [];
par.phi      = [];
par.psi      = [];
par.theta    = [];
par.vpsi     = [];
par.vphi     = [];
par.vtheta   = [];
par.array{1} = [];
par.arraydim = [];
par.seqcon   = [];
par.lro      = [];
par.lpe      = [];
par.lpe2     = [];
par.pro      = [];
par.ppe      = [];
par.ppe2     = [];
par.gap      = [];
par.pss      = [];
par.thk      = [];
par.thk2     = [];
par.pos1     = [];
par.pos2     = [];
par.pos3     = [];
par.vox1     = [];
par.vox2     = [];
par.vox3     = [];
par.nD       = [];
par.cntr     = [];
par.gain     = [];
par.shimset  = [];
par.z1       = [];
par.z2       = [];
par.z3       = [];
par.z4       = [];
par.z5       = [];
par.z6       = [];
par.z7       = [];
par.z8       = [];
par.x1       = [];
par.y1       = [];
par.xz       = [];
par.yz       = [];
par.xy       = [];
par.x3       = [];
par.y3       = [];
par.x4       = [];
par.y4       = [];
par.z1c      = [];
par.z2c      = [];
par.z3c      = [];
par.z4c      = [];
par.xz2      = [];
par.yz2      = [];
par.xz2      = [];
par.yz2      = [];
par.zxy      = [];
par.z3x      = [];
par.z3y      = [];
par.zx3      = [];
par.zy3      = [];
par.z4x        = [];
par.z4y        = [];
par.z5x        = [];
par.z5y        = [];
par.x2y2       = [];
par.z2xy       = [];
par.z3xy       = [];
par.z2x3       = [];
par.z2y3       = [];
par.z3x3       = [];
par.z3y3       = [];
par.z4xy       = [];
par.zx2y2      = [];
par.z2x2y2     = [];
par.z3x2y2     = [];
par.z4x2y2     = [];
par.petable    = [];
par.petype     = [];
par.eta        = [];
par.nrcvrs     = [];
par.trise      = [];
par.at         = [];
par.gro        = [];
par.gmax       = [];
par.intlv      = [];
par.rcvrs      = [];
par.celem      = [];
par.arrayelemts = [];
par.contrast   = [];
par.eff_echo   = [];
par.tep        = [];
par.date       = [];
par.ti         = [];
par.ti2        = [];
par.gss2       = [];
par.gss        = [];
par.tpwri      = [];
par.tpwr1      = [];
par.tpwr2      = [];
par.tpwrf1     = [];
par.tpwrf2     = [];
par.orient     = [];
par.tof        = [];
par.resto      = [];
par.grox       = [];
par.groy       = [];
par.fov        = [];
par.res        = [];
par.npix       = [];
par.nseg       = [];
par.nzseg      = [];
par.waveform   = [];
par.SR         = [];
par.gradfrac   = [];
par.sfrq       = [];
par.B0         = [];
par.dtmap      = [];
par.nnav       = [];
par.tnav       = [];
par.fast       = [];
par.bt         = [];
par.nhomo      = [];
par.fpmult     = [];
par.d1         = [];
par.ss         = [];
par.ssc        = [];
par.r1         = [];
par.r2         = [];
par.ps_coils   = [];
par.coil_array = [];
par.nav        = [];
par.fliplist   = [];
par.varflip    = [];
par.fliptrain  = [];
par.spiral_gmax= [];
par.nfreq      = [];
par.rec_flag   = [];
par.flip       = [];
par.flipprep   = [];
par.seg        = [];
par.state      = [];
par.rfdelay    = [];
par.gro        = [];
par.gimp       = [];
par.SR         = [];
par.readaxis   = [];
par.etl        = [];
par.grof       = [];
par.grof2      = [];
par.Po         = [];
par.Psl        = [];
par.comment    = [];
par.haste      = [];
par.hasteLines = [];
par.name       = [];
par.pks        = [];
par.fcorr      = [];
par.fcorrx     = [];
par.fcorrt1    = [];
par.fcorrt2    = [];
par.ffilt      = [];
par.fftype     = [];
par.ffcut      = [];
par.interp     = [];
par.interpnp   = [];
par.interpnv   = [];
par.interppss  = [];
par.dcrmv      = [];
par.tab_converted = [];
par.flash_converted=[];
par.drec       = [];
par.senseR     = [];
par.sloop      = [];
par.ir         = [];
par.rfcoil     = [];
par.time_complete=[];
par.pocs       = [];
par.arrayrec   = [];
par.pss0       = [];
par.t2map      = [];
par.tep        = [];
par.petype     = [];
par.nt2        = [];
par.nav_enable = [];
par.gpat       = [];
par.shots      = [];
par.spinecho   = [];
par.sp_in      = [];
par.noGref     = [];
par.B0ref      = [];
par.nt2a       = [];
par.ref        = [];
par.gcrush     = [];
par.tcrush     = [];
par.tspoil     = [];
par.gspoil     = [];
par.np_a       = [];
par.nv_a       = [];
par.nv_c       = [];
par.NS         = [];
par.RO_CVRD    = [];
par.studyid_   = [];

% The structure par MUST use parameter names which match the
% parameter names used by VNMR or it will not read correctly

% This function extracts the first letter of each field name
% in order to quickly discard any line which does not being 
% with one of the parameters of interest.
names = fieldnames(par);
value = names{1}(1);
for n = 2:size(names,1)
   if isempty(findstr(value, names{n}(1)))
      value = [value names{n}(1)];
   end
end
clear n names;

buffer = fgets(fp);
 
% Parse the ASCII procpar file
while ( buffer ~= -1 )
   % Check to see if the first letter in the buffer matches
   % the first character of any parameter of interest
   % This provides a 3 to 4 fold speed up.
   if (findstr(value, buffer(1)))
      
      % Get only the first word of the buffer, the parameter name
      ind = findstr(buffer, ' ');
      lenb = ind(1)-1;
      buffer = buffer(1:lenb);
      
      % Read in required parameters
      if (lenb == 2)
         if any(strcmp(buffer, {'z1','z2','z3','z4','z5','z6','z7','z8', ...
               'x1','y1','xz','yz','xy','x3','y3','x4','y4','nD', ...
               'nf','ni','np','ns','nv','ne','ss','nt','ti','te', ...
               'tr','sw','at','bt','d1','SR','r1','r2','B0','SR', ...
               'Po','NS'}))
            val = sscanf(fgets(fp), '%f');
            par.(buffer) = val(2:end).';
            fgets(fp);
         elseif any(strcmp(buffer, {'dp','ir'}))
            tmpbuffer = fgets(fp);
            ind = findstr(tmpbuffer, '"');
            par.(buffer) = tmpbuffer(ind(1)+1:ind(2)-1);
            fgets(fp);
         end
      elseif (lenb == 3),
         if any(strcmp(buffer, {'z1c','z2c','z3c','z4c','xz2','yz2', ...
               'zxy','z3x','z3y','zx3','zy3','z4x','z4y','z5x','z5y', ...
               'ssc','nv2','nv3','ne2','nss','nse','tep','esp','lro', ...
               'pro','lpe','ppe','gap','pss','phi','psi','thk','gro', ...
               'gss','tof','fov','res','gro','etl','Psl','ti2','eta', ...
               'nt2','tep','ref'}))
            val = sscanf(fgets(fp), '%f');
            par.(buffer) = val(2:end).';
            fgets(fp);
         elseif any(strcmp(buffer, {'nav','seg','pks'}))
            tmpbuffer = fgets(fp);
            ind = findstr(tmpbuffer, '"');
            par.(buffer) = tmpbuffer(ind(1)+1:ind(2)-1);
            fgets(fp);
         end
      elseif (lenb == 4)
         if any(strcmp(buffer, {'x2y2','z2xy','z3xy','z2x3','z2y3', ...
               'z3x3','z3y3','z4xy','lpe2','ppe2','thk2','pos1', ...
               'pos2','pos3','thk2','vox1','vox2','vox3','vpsi', ...
               'vphi','gmax','flip','gss2','sfrq','nnav','tnav', ...
               'grox','groy','nseg','npix','gain','cntr','gimp', ...
               'grof','drec','pss0','nt2a','nv_a','np_a','nv_c'}))
            val = sscanf(fgets(fp), '%f');
            par.(buffer) = val(2:end).';
            fgets(fp);
         elseif any(strcmp(buffer,{'date','fast','name','pocs','gpat'}))
            tmpbuffer = fgets(fp);
            ind = findstr(tmpbuffer, '"');
            par.(buffer) = tmpbuffer(ind(1)+1:ind(2)-1);
            fgets(fp);
         end
      elseif (lenb == 5)
         if any(strcmp(buffer, {'zx2y2','celem','tproj','trise', ...
               'tpwr1','tpwr2','tpwri','theta','resto','nhomo', ...
               'nfreq','dtmap','nzseg','nproj','state','ffcut', ...
               't2map','grof2','shots','sp_in','B0ref'}))
            val = sscanf(fgets(fp), '%f');
            par.(buffer) = val(2:end).';
            fgets(fp);
         elseif any(strcmp(buffer,{'haste','fcorr','ffilt','dcrmv','sloop'}))
            tmpbuffer = fgets(fp);
            ind = findstr(tmpbuffer, '"');
            par.(buffer) = tmpbuffer(ind(1)+1:ind(2)-1);
            fgets(fp);
         elseif (strcmp(buffer, 'rcvrs'))
            buffer = fgets(fp);
            % Count the number of 'y' values to find number of receivers
            ind = findstr(buffer, '"');
            buffer = buffer(ind(1)+1:ind(2)-1);
            par.rcvrs = buffer;
            par.nrcvrs = sum(double(buffer) == 121);
            fgets(fp);
         elseif (strcmp(buffer, 'intlv'))
            buffer = fgets(fp);
            ind = findstr(buffer, '"');
            par.intlv = buffer(ind(1)+1:ind(2)-1);
            fgets(fp);
         elseif (strcmp(buffer, 'array'))
            buffer = fgets(fp);
            % Strip the buffer down to the core data
            ind = findstr(buffer, '"');
            buffer = buffer((ind(1)+1):(ind(2)-1));
            
            % Parse the data string
            index1 = 1;
            index2 = 1;
            incr = 0;
            par.array{1}{1} = [];
            for o = 1:length(buffer)
               switch buffer(o)
               case '('
                  incr = incr + 1;
                  if incr > 1
                     error('invalid coupling in array paramater.');
                  end
               case ')'
                  incr = incr - 1;
                  if incr == 0
                     index2 = 1;
                  elseif incr < 0
                     error('invalid coupling in array paramater.');
                  end
               case ','
                  if incr > 0
                     index2 = index2 + 1;
                  else
                     index1 = index1 + 1;
                  end
                  par.array{index1}{index2} = [];
               otherwise
                  par.array{index1}{index2} = [par.array{index1}{index2} buffer(o)];
               end
            end
            fgets(fp);
         end
      elseif (lenb == 6)
         if any(strcmp(buffer, {'z2x2y2','z3x2y2','z4x2y2','vtheta', ...
               'fpmult','senseR','tpwrf1','tpwrf2','noGref','gcrush',...
               'tcrush','tspoil','gspoil'}))
            val = sscanf(fgets(fp), '%f');
            par.(buffer) = val(2:end).';
            fgets(fp);
         elseif any(strcmp(buffer, {'seqfil','seqcon','orient','fcorrx','fftype','interp','rfcoil','petype'})) 
            tmpbuffer = fgets(fp);
            ind = findstr(tmpbuffer, '"');
            par.(buffer) = tmpbuffer(ind(1)+1:ind(2)-1);
            fgets(fp);
         end
      elseif (lenb == 7)
         if any(strcmp(buffer, {'espincr','rfphase','shimset','rfdelay','fcorrt1','fcorrt2','RO_CVRD'}))
            val = sscanf(fgets(fp), '%f');
            par.(buffer) = val(2:end).';
            fgets(fp);
         elseif any(strcmp(buffer, {'rfspoil','varflip','petable','comment'}))
            tmpbuffer= fgets(fp);
            ind = findstr(tmpbuffer, '"');
            par.(buffer) = tmpbuffer(ind(1)+1:ind(2)-1);
            fgets(fp);
         end
      elseif (lenb == 8)
         if any(strcmp(buffer, {'evenecho','polarity','arraydim', ...
               'fliplist','gradfrac','eff_echo','interpnp','interpnv'}))
            val = sscanf(fgets(fp), '%f');
            par.(buffer) = val(2:end).';
            fgets(fp);
         elseif any(strcmp(buffer, {'waveform','contrast', ...
               'flipprep','readaxis','rec_flag','arrayrec','spinecho',...,
               'studyid_'}))
            tmpbuffer = fgets(fp);
            ind = findstr(tmpbuffer, '"');
            par.(buffer) = tmpbuffer(ind(1)+1:ind(2)-1);
            fgets(fp);
         elseif (strcmp(buffer, 'ps_coils'))
            tmpbuffer = fgets(fp);
            n = sscanf(tmpbuffer, '%f', 1);
            for m = 1:n
               ind = findstr(tmpbuffer, '"');
               par.(buffer){m} = tmpbuffer(ind(1)+1:ind(2)-1);
               tmpbuffer = fgets(fp);
            end
         end
      elseif (lenb == 9)
          if any(strcmp(buffer,{'fliptrain','interppss'}))
              val = sscanf(fgets(fp),'%f');
              par.(buffer) = val(2:end).';
              fgets(fp);
          end
      elseif (lenb == 10)
         if any(strcmp(buffer, {'coil_array','nav_enable'}))
            tmpbuffer = fgets(fp);
            n = sscanf(tmpbuffer, '%f', 1);
            for m = 1:n
               ind = findstr(tmpbuffer, '"');
               par.(buffer){m} = tmpbuffer(ind(1)+1:ind(2)-1);
               tmpbuffer = fgets(fp);
            end
         elseif any(strcmp(buffer,{'hasteLines'}))
             val = sscanf(fgets(fp),'%f');
             par.(buffer) = val(2:end).';
             fgets(fp);
         end
      elseif (lenb == 11)
         if any(strcmp(buffer, {'arrayelemts','spiral_gmax'}))
            val = sscanf(fgets(fp), '%f');
            par.(buffer) = val(2:end).';
            fgets(fp);
         end
      elseif (lenb == 13)
          if any(strcmp(buffer,{'tab_converted'}))
              val = sscanf(fgets(fp),'%f');
              par.(buffer) = val(2:end).';
              fgets(fp);
          elseif any(strcmp(buffer, {'time_complete'}))
            tmpbuffer = fgets(fp);
            ind = findstr(tmpbuffer, '"');
            par.(buffer) = tmpbuffer(ind(1)+1:ind(2)-1);
            fgets(fp);
          end
      elseif (lenb == 15)
          if any(strcmp(buffer,{'flash_converted'}))
              val = sscanf(fgets(fp),'%f');
              par.(buffer) = val(2:end).';
              fgets(fp);
          end
      end
   end

   buffer = fgets(fp);
end  % End of while loop

% Close file
fclose(fp);

return
