datax=get(hedittest1,'String'); 
%eval([mCallback]);
%figure(fig)
%axis('off');
hti=1;
%eval(['[Dn Dnp Dnm p test]=kuiper(' datax ');'])
eval(['[Vn K V]=kuiper(' datax ');'])
test='NS';
if V>1.747, test='*'; end
if V>2.001, test='**'; end
ht(hti)=text(0.02,0.44,'Dn:'); hti=hti+1;
hedit=uicontrol('Style','edit','String',num2str(Vn),'Position',POSITION(1,:));
set(hedit,'BackgroundColor',[1 1 1]);
ht(hti)=text(0.02,0.32,'Dnm:'); hti=hti+1;
hedit=uicontrol('Style','edit','String',num2str(K),'Position',POSITION(2,:));
set(hedit,'BackgroundColor',[1 1 1]);
ht(hti)=text(0.02,0.2,'Dnp:'); hti=hti+1;
hedit=uicontrol('Style','edit','String',num2str(V),'Position',POSITION(3,:));
set(hedit,'BackgroundColor',[1 1 1]);
%ht(hti)=text(0.53,0.32,'P:'); hti=hti+1;
%hedit=uicontrol('Style','edit','String',num2str(p),'Position',[300 85 80 20]);
%set(hedit,'BackgroundColor',[1 1 1]);
ht(hti)=text(0.53,0.2,'test:'); hti=hti+1;
hedit=uicontrol('Style','edit','String',test,'Position',[300 60 80 20]);
set(hedit,'BackgroundColor',[1 1 1]);
%set(ht,'FontName','Serif');
%set(ht,'FontSize',27);
%%uicontrol('Style','pushbutton','String','Plot','Callback','sg6p4p1b','Position',[360 200 60 27]);
%uicontrol('Style','pushbutton','String','Table','Callback','sg6p4p1c','Position',[360 170 60 27]);

