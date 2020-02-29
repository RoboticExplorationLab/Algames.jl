export
    test_allocation

function test_allocation(f, args)
    @allocated f(args...)
end

#
#
# @test true
# # Test Passed
#
# @test [1, 2] + [2, 1] == [3, 3]
# # Test Passed
#
# @test_throws BoundsError [1, 2, 3][4]
# # Test Passed
#       # Thrown: BoundsError
#
# @test_throws DimensionMismatch [1, 2, 3] + [1, 2]
# # Test Passed
#       # Thrown: DimensionMismatch
