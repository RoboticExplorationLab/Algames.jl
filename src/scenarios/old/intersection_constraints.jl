# Add the intersection constraints to the car with id player_id driving on lane ∈ [1,4].
function add_intersection_constraints(conSet::ConstraintSet, l1, l2, circle_radius, car_radius,
    player_id, px, N; lane=1)
    P = landmarks(car_radius, l1-car_radius, l2)
    r = circle_radius + car_radius
    if lane == 1
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[1],  P[6],  P[5],  r,
            player_id, px, N)
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[2],  P[19], P[22], r,
            player_id, px, N)
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[28], P[25], P[30], r,
            player_id, px, N)
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[29], P[12], P[11],  r,
            player_id, px, N)
    elseif lane == 2
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[3],  P[8],  P[5],  r,
            player_id, px, N)
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[4],  P[9],  P[10], r,
            player_id, px, N)
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[16], P[13], P[30], r,
            player_id, px, N)
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[29], P[12], P[11], r,
            player_id, px, N)
    elseif lane == 3
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[3],  P[8],  P[5],  r,
            player_id, px, N)
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[22], P[21], P[4],  r,
            player_id, px, N)
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[28], P[27], P[32], r,
            player_id, px, N)
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[31], P[14], P[11],  r,
            player_id, px, N)
    elseif lane == 4
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[3],  P[20], P[17],  r,
            player_id, px, N)
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[4],  P[21], P[22], r,
            player_id, px, N)
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[28], P[25], P[30], r,
            player_id, px, N)
        add_rounded_boundary_constraints(conSet::ConstraintSet,
            P[29], P[24], P[23],  r,
            player_id, px, N)
    end
    return nothing
end

function landmarks(l1, l2, l3)
    P = [[-l2, l3], #Line1
         [-l1, l3],
         [ l1, l3],
         [ l2, l3],
         [-l3, l2], #Line2
         [-l2, l2],
         [-l1, l2],
         [ l1, l2],
         [ l2, l2],
         [ l3, l2],
         [-l3, l1], #Line3
         [-l2, l1],
         [-l1, l1],
         [ l1, l1],
         [ l2, l1],
         [ l3, l1],
         [-l3,-l1], #Line4
         [-l2,-l1],
         [-l1,-l1],
         [ l1,-l1],
         [ l2,-l1],
         [ l3,-l1],
         [-l3,-l2], #Line5
         [-l2,-l2],
         [-l1,-l2],
         [ l1,-l2],
         [ l2,-l2],
         [ l3,-l2],
         [-l2,-l3], #Line6
         [-l1,-l3],
         [ l1,-l3],
         [ l2,-l3]]
    return P
end
