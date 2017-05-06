typedef struct
{
  double  ***x;
  double  ***y;
  double  ***z;
} Vector;

typedef struct
{
  double  ***xx;
  double  ***yy;
  double  ***zz;
  double  ***xy;
  double  ***xz;
  double  ***yz;
} Tensor;

typedef struct
{
  double ***s;
} Scalar;

typedef struct
{
  double ***xxx;
  double ***xxy;
  double ***xxz;
  double ***xyy;
  double ***xyz;
  double ***xzz;
  double ***yxx;
  double ***yxy;
  double ***yxz;
  double ***yyy;
  double ***yyz;
  double ***yzz;
  double ***zxx;
  double ***zxy;
  double ***zxz;
  double ***zyy;
  double ***zyz;
  double ***zzz;
} Connection;

