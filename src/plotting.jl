using .MATLAB

export trisurf, double_trisurf, trisurf_curvature

function trisurf(p, t, fignum = 0; edgecolor = "k", linewidth = 0.5, facecolor = [])    
    if fignum == 0
        mat"figure()"
    else
        mat"figure(double($fignum))"
    end
    
    if isempty(facecolor)
        mat"trisurf($t', $p(:,1), $p(:,2), $p(:,3), 'edgecolor', 'none')"
        mat"shading interp; hold on;"
        mat"trisurf($t', $p(:,1), $p(:,2), $p(:,3), 'edgecolor', $edgecolor, 'LineWidth', $linewidth, 'facecolor', 'none')"
        mat"hold off;"
    else
        mat"trisurf($t', $p(:,1), $p(:,2), $p(:,3), 'edgecolor', $edgecolor, 'LineWidth', $linewidth, 'facecolor', $facecolor)"
    end
    
    mat"axis equal;"
    mat"camlight('right')"
    mat"material([0.4 0.8 0])"
end

function trisurf(p::Vector, t, fignum = 0; edgecolor = "k", linewidth = 0.5, facecolor = [])    
    tmp = permutedims(reshape(p, 3, length(p) ÷ 3), (2,1))
    
    trisurf(tmp, t, fignum; edgecolor = edgecolor, linewidth = linewidth, facecolor = facecolor)
end

function double_trisurf(x1::Vector, x2::Vector, t1, t2, fignum = 0)    
    if fignum == 0
        mat"figure()"
    else
        mat"figure(double($fignum))"
    end
    
    tmp1 = permutedims(reshape(x1, 3, length(x1) ÷ 3), (2,1))
    tmp2 = permutedims(reshape(x2, 3, length(x2) ÷ 3), (2,1))

    mat"trisurf($t1', $tmp1(:,1), $tmp1(:,2), $tmp2(:,3), 'edgecolor', 'r', 'facecolor', 'r')"
    mat"hold on;"
    mat"trisurf($t2', $tmp2(:,1), $tmp2(:,2), $tmp2(:,3), 'edgecolor', 'k', 'facecolor', 'k')"
    mat"hold off;"
    
    mat"axis equal;"
end

function trisurf_curvature(p::Vector, mi::MeshInfo, fignum = 0)    
    if fignum == 0
        mat"figure()"
    else
        mat"figure(double($fignum))"
    end
    
    t = mi.t
    H,K,N,A = calculate_curvature(p, mi::MeshInfo)
    
    mat"trisurf($t', $p(1:3:end), $p(2:3:end), $p(3:3:end), $H, 'edgecolor', 'none')"
    mat"shading interp; hold on;"
    mat"hold off;"
    mat"colorbar"
    
    mat"axis equal;"
    mat"camlight('right')"
    mat"material([0.4 0.8 0])"
end

function add_text_labels(p::Vector)
    mat"hold on"

    for i = 1:length(p)÷3
        mat"text($p(3*$i-2), $p(3*$i-1), $p(3*$i)+0.05, sprintf('%d', $i))"
    end


    mat"hold off"
end

function plot_octree(o::MeshOctNode)
    mat"hold on"

    q = [o]

    while !isempty(q)
        node = popfirst!(q)
        append!(q, children(node))

        node_rect = Rect(node)
        (rx,ry,rz) = node_rect.x
        (wx,wy,wz) = node_rect.w

        mat"plot3([$rx, $rx+$wx, $rx+$wx, $rx, $rx], [$ry, $ry, $ry+$wy, $ry+$wy, $ry], [$rz, $rz, $rz, $rz, $rz], 'k-')"
        mat"plot3([$rx, $rx+$wx, $rx+$wx, $rx, $rx], [$ry, $ry, $ry+$wy, $ry+$wy, $ry], [$rz, $rz, $rz, $rz, $rz]+$wz, 'k-')"
        mat"plot3([$rx, $rx], [$ry, $ry], [$rz, $rz+$wz], 'k-')"
        mat"plot3([$rx+$wx, $rx+$wx], [$ry, $ry], [$rz, $rz+$wz], 'k-')"
        mat"plot3([$rx+$wx, $rx+$wx], [$ry+$wy, $ry+$wy], [$rz, $rz+$wz], 'k-')"
        mat"plot3([$rx, $rx], [$ry+$wy, $ry+$wy], [$rz, $rz+$wz], 'k-')"
    end

    mat"hold off"
end
