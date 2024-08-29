% Unit tests for the preprocessing class
classdef test_preprocessing < matlab.unittest.TestCase

    methods(Test)

        % Test preprocessing.set_nfcts
        function test_set_nfcts(testCase)
            expmnt = preprocessing(102, '../examples/');
            preprocessing.set_nfcts(expmnt, 666);
            testCase.verifyEqual(expmnt.nfcts, 666);
        end
    end
end
