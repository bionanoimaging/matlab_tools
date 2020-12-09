% aViewer=LiveView(aViewer) : updates a live viewer.
% 
function aViewer=LiveView(aViewer,aRecon)
if isempty(aViewer) || (isnumeric(aViewer) && aViewer > 0)
    dipshow(3,aRecon);drawnow();
else
    if aViewer == 0
        aViewer = view5d(aRecon); 
    else
        try
            aViewer.toFront();
            view5d(aRecon,0,'replaceElement',aViewer,0);
            aViewer.ProcessKeyMainWindow('t');
            aViewer.UpdatePanels()
        catch
            dipshow(3,aRecon);
            drawnow();
        end
    end
end
