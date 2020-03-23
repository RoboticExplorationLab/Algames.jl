export
    test_allocation

function test_allocation(f, args)
    @allocated f(args...)
end
