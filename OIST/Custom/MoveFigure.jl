#!/bin/julia
function MoveFigure(f,x,y)
#   Move figure's upper left corner to pixel (x, y)
    backend = matplotlib.get_backend()
    if backend == "TkAgg"
        newp = "+$x+$y"
        f.canvas.manager.window.wm_geometry(newp)
    elseif backend == "WXAgg"
        f.canvas.manager.window.SetPosition((x, y))
    else
        # This works for QT and GTK
        # You can also use window.setGeometry
        f.canvas.manager.window.move(x, y)
    end
end
#----------------------------------------------------------------------   
