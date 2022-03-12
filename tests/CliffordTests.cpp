#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "Clifford.h"

TEST_SUITE("Clifford Class Tests")
{

    TEST_CASE("p and q are assigned correctly.")
    {
        Clifford cliff11(1, 1);
        CHECK(cliff11.get_q() == 1);
        CHECK(cliff11.get_q() == 1);
    }

    TEST_CASE("matrix dim is set correctly")
    {
        SUBCASE("1-dim Clifford algebras are correct size.")
        {
            SUBCASE("Type (1,0) is 1x1")
            {
                // 1-dimensionals
                Clifford cliff10(1, 0);
                CHECK(cliff10.get_matrix_size() == 1);
            }
            SUBCASE("Type (0,1) is 1x1")
            {
                Clifford cliff01(0, 1);
                CHECK(cliff01.get_matrix_size() == 1);
            }
        }

        SUBCASE("2-dim Clifford are correct size.")
        {
            // 2-dimensionals
            SUBCASE("Type (1,1) is 2x2 ")
            {
                Clifford cliff11(1, 1);
                CHECK(cliff11.get_matrix_size() == 2);
            }
            SUBCASE("Type (0,2) is 2x2")
            {
                Clifford cliff02(0, 2);
                CHECK(cliff02.get_matrix_size() == 2);
            }

            SUBCASE("Type (2,0) is 2x2")
            {
                Clifford cliff20(2, 0);
                CHECK(cliff20.get_matrix_size() == 2);
            }
        }

        SUBCASE("3-dim Clifford algebras are correct size")
        {
            SUBCASE("Type(3,0) is 2x2")
            {
                // 3-dimensonals
                Clifford cliff30(3, 0);
                CHECK(cliff30.get_matrix_size() == 2);
            }

            SUBCASE("Type(0,3) is 2x2")
            {
                Clifford cliff03(0, 3);
                CHECK(cliff03.get_matrix_size() == 2);
            }

            SUBCASE("Type(1,2) is 2x2")
            {
                Clifford cliff12(1, 2);
                CHECK(cliff12.get_matrix_size() == 2);
            }
            SUBCASE("Type (2,1) is 2x2")
            {
                Clifford cliff21(2, 1);
                CHECK(cliff21.get_matrix_size() == 2);
            }
        }

        // 4-dimensionals
        SUBCASE("Type (1,3) is 4x4")
        {
            Clifford cliff13(1, 3);
            CHECK(cliff13.get_matrix_size() == 4);
        }
    }
}