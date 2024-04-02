function [IE,IEm,TE] = IE_filter(pts, twindow)

% FIFO
IE_k = zeros(1280,720+2); % positive fifo - index
Pt   = zeros(1280,720+2); % positive fifo - index
Pp   = zeros(1280,720+2); % positive fifo - index


%  initialize
IE = false(length(pts.x),1);
TE = false(length(pts.x),1);
IEm = zeros(length(pts.x),1);

for k=1:length(pts.x)
        
        
        if (pts.p(k)~=Pp(pts.x(k),pts.y(k))) || (pts.ts(k)-Pt(pts.x(k),pts.y(k))>twindow) 
            % (polarity switch) or (long time since last event) => new inceptive event
            
            IE(k)=true; % mark current
            IEm(k)=1; % magnitude of 1
            IE_k(pts.x(k),pts.y(k)) = k; % remember new IE location
            
        else % trailing event

            TE(k) = true; % mark current
            IEm(IE_k(pts.x(k),pts.y(k)))=IEm(IE_k(pts.x(k),pts.y(k)))+1; % increment IE magnitude 

        end

        % update
        Pt(pts.x(k),pts.y(k))=pts.ts(k); % previous time stamp
        Pp(pts.x(k),pts.y(k))=pts.p(k); % previous polarity

end