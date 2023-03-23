class Test_Gyro : public ::testing::Test
{
protected:
    void SetUp() override
    {
        gyro::init_params();
    }
}

TEST_F(Test_Gyro, distBoundary)
{
    gyro::dcntl[Param::R0] = 0.1;
    dcntl[Param::R]        = 1;
    EXPECT_EQ(gyro::distBoundary(0.5, 1, 1), 0.2);
    EXPECT_EQ(gyro::distBoundary(0.5, 0.342, 3), 0.2);
    EXPECT_EQ(gyro::distBoundary(0.1, 1, 1), 0);
    EXPECT_EQ(gyro::distBoundary(1, 1, 1), 0);

    //To ask: What if DirBC_Interior = 0 ? If we discretize across the origin then the interior circle
    //is not part of the dirichlet boundary is it not ?
}

TEST(Test_Gyro, sign)
{
    EXPECT_EQ(gyro::sign(1), 1);
    EXPECT_EQ(gyro::sign(0), 1);
    EXPECT_EQ(gyro::sign(-45), -1);
    EXPECT_EQ(gyro::sign(642), 1);
}

/*TEST(Test_Gyro, select_functions_class)
{
    gyro::select_functions_class(0,0,0,5); //circular geometry
    gyro::functions ->
}


TEST_F(Test_Gyro, eval_sol){
    gyro::select_functions_class(0,0,0,5); //polar coord
    gyro::functions -> 
}*/