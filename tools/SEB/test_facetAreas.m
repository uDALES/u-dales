% Unit tests for facetAreas.m
classdef test_facetAreas < matlab.unittest.TestCase

    methods(Test, TestTags = {'Unit'})

        % Test facetAreas
        function test(testCase)
            F = [
                    [1,2,3];
                    [3,1,2];
                    [2,3,1]
                ];
            V = magic(3);
            actual = facetAreas(F, V);
            expected = [20.7846; 20.7846; 20.7846];

            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-4);
        end
    end
end