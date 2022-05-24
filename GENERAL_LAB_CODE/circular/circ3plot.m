function circplot(M,z,deg,nr,tot,circ)
%
% Plots directional data in a circular diagram.
%
% CALL: circplot(x,z,deg,nr,tot)
%
% where
%
%     x   = column vector of directions,
%     deg = 'deg' if measured in degrees.
%     nr  = (optional) position of the plot (default=1 plot)
%     tot = (optional) total number of plots (default=1 or 6 (if nr>0))
%
% Copyright (C) 1993, Bjorn Holmquist, Dept of Math. Stat., University of Lund.
%
global htext
global htextcount
global FontSize
if nargin<5, circ=1; end
if (nargin==3),
 if (nr==0), subplot(1,1,1)
 else
  subplot(2,3,nr)
 end
end
if (nargin==4),
 if (nr==0), subplot(1,1,1)
 else
  if (tot>=9),
    subplot(3,3,nr)
  end
  if ((tot>6)&(tot<9)),
    subplot(2,4,nr)
  end
  if ((tot>=5)&(tot<7)),
    subplot(2,3,nr)
  end
  if ((tot>=3)&(tot<5)),
    subplot(2,2,nr)
  end
  if (tot==2),
    subplot(1,2,nr)
  end
 end
end
if (nargin==2), nr=0; end
if (nargin>0)
  uppl=1500;
  X=cos(2*pi*[0:uppl]/uppl);
  Y=sin(2*pi*[0:uppl]/uppl);
  hold on
%  axis([-1.5,1.5,-1.5,1.5]);
  axis([-1.6,1.6,-1.6,1.6]);
  axis('off')
  axis('equal') 
  if circ~=0, plot(X,Y,'k'); end 
  if (deg(1:3)=='rad')
     [sortord,y]=sort(M(:,1));
		 z=z(sortord);
     [ry,cy]=size(y);
     r1=ones(ry,1);
     for i=2:ry,
       if (y(i)-y(i-1)<0.03),
         r1(i,1)=r1(i,1)+0.03;
       end
     end
     plot3(r1.*cos(pi/2-y),r1.*sin(pi/2-y),z,'bo');
  else
    qua=[0,90,180,270]; co=cos(pi*(90-qua)/180); si=sin(pi*(90-qua)/180);
    plot([1.0*co;0.95*co],[1.0*si;0.95*si],'k-')
  end
  y=sort(M(:,1));
  [ry,cy]=size(y);
  r1=ones(ry,1)+0.05;
  for i=2:ry,
    ant=sum(abs(y(i)*ones(i-1,1)-y([1:i-1]))<1);
    r1(i,1)=r1(i,1)+ant*0.05;
  end
end
  plot3(r1.*cos(pi*(90-y)/180),r1.*sin(pi*(90-y)/180),z,'bo')
%   htextcount=htextcount+1;
%   htext(htextcount)=text(0,1.1,'^','HorizontalAlignment','center'); 
  htextcount=htextcount+1;
  htext(htextcount)=text(0,1.2,'N','HorizontalAlignment','center'); %,'Rotation',90)
  htextcount=htextcount+1;
  htext(htextcount)=text(1.2,0,'E','VerticalAlignment','middle','HorizontalAlignment','center');
  htextcount=htextcount+1;
  htext(htextcount)=text(0,-1.2,'S','HorizontalAlignment','center');
  htextcount=htextcount+1;
  htext(htextcount)=text(-1.2,0,'W','VerticalAlignment','middle','HorizontalAlignment','center');
  %if (nr>0), 
  set(htext,'FontSize',FontSize); %end
  hold on
%end


