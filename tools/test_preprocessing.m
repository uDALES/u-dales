% Unit tests for the preprocessing class
classdef test_preprocessing < matlab.unittest.TestCase

    methods(Test, TestTags = {'Unit'})

        % Test preprocessing.addvar
        function test_addvar(testCase)
            expmnt = preprocessing(102, '../examples/');

            % Add a new variable via addvar, check its value
            preprocessing.addvar(expmnt, "a_new_variable", 123);
            testCase.verifyEqual(expmnt.a_new_variable, 123);

            % Try adding the same variable via addvar, check its value
            preprocessing.addvar(expmnt, "a_new_variable", 456);
            testCase.verifyNotEqual(expmnt.a_new_variable, 456);
        end

        % Test preprocessing.set_nfcts
        function test_set_nfcts(testCase)
            expmnt = preprocessing(102, '../examples/');
            preprocessing.set_nfcts(expmnt, 666);
            testCase.verifyEqual(expmnt.nfcts, 666);
        end
    end
end
