%[R,n,p,h,test]=kruswall(y,xc);
[l,p,test,df]=bartlett(y,xc)

fbartlett=figure('Name','Bartlett''s test','NumberTitle','off','Position',[12 12 390 250],'Resize','off');
set(f6p3p1,'MenuBar','none');axis('off'); 
yip=[185,155,125,90,60,30];
uicontrol('Style','frame','Position',[18 15 295 206],'Backgroundcolor',Bgc);
uicontrol('Style','pushbutton','String','OK','Callback','close(fbartlett)','Position',[340 15 40 27],'BackgroundColor',OKBgc);


%uicontrol('Style','text','String',num2str(R),'Position',[75 yip(3) 180 20],'BackgroundColor',Bgc2);
uicontrol('Style','text','String',num2str(l),'Position',[175 yip(4) 80 20],'BackgroundColor',Bgc2);
uicontrol('Style','text','String',num2str(p),'Position',[175 yip(5) 80 20],'BackgroundColor',Bgc2);
uicontrol('Style','text','String',test,'Position',[175 yip(6) 80 20],'BackgroundColor',Bgc2);
%uicontrol('Style','text','String','Lambda'':','Position',[20 yip(3) 60 20],'HorizontalAlignment','Left','BackgroundColor',Bgc,'ForegroundColor',Fgc);
uicontrol('Style','text','String','Lambda'':','Position',[20 yip(4) 150 20],'HorizontalAlignment','Left','BackgroundColor',Bgc,'ForegroundColor',Fgc);
uicontrol('Style','text','String','Test:','Position',[20 yip(6) 150 20],'HorizontalAlignment','Left','BackgroundColor',Bgc,'ForegroundColor',Fgc);
uicontrol('Style','text','String',['Probability P(Lambda''>chi2(' num2str(df) ')):'],'Position',[20 yip(5) 155 20],'HorizontalAlignment','Left','BackgroundColor',Bgc,'ForegroundColor',Fgc);
