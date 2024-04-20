 double CalculateFractalDimension(int partID,double porosity,int genus)
{
    std::vector<my::render::ViewData> viewdatas;
    model::AppImp::Instance()->_model->getPartById(partID)->getViewData(viewdatas);
    int i, j;
    std::vector<double> pntCoords;
    std::vector<int> triangle3Vers;
    for (i = 0; i < viewdatas.size(); ++i)
    {
//        if (!viewdatas[i].mesh)
//            continue;
        if (viewdatas[i].trias.empty())
            continue;
        triangle3Vers.insert(triangle3Vers.end(), viewdatas[i].trias.begin(), viewdatas[i].trias.end());
        for (j = 0; j < viewdatas[i].nodes.size(); ++j)
            pntCoords.push_back(viewdatas[i].nodes[j]);
    }
    if (triangle3Vers.empty())
    {
        Messager::Msg(Messager::eWarning, "Invalid part or no mesh data in this part!");
        return -1.0;
    }

    mTranslate obj;
    double fractalDimension = obj.CalcFractalDimension(pntCoords.data(), pntCoords.size() / 3, triangle3Vers.data(), triangle3Vers.size() / 3, porosity, genus, true);
    if(fractalDimension<0.0)
    {
        Messager::Msg(Messager::eWarning, "Can't calculate fractal dimension!");
        return -1.0;
    }

    double minVal = 2.0;
    double maxVal = 3.0;
    double numOfFD = fractalDimension;

    fractalDimension = (fractalDimension - minVal) / (maxVal - minVal) * 100.0;

    return fractalDimension;
}

extern long long FillTriangle(const double ver1[3], const double ver2[3], const double ver3[3],
    const double boxMinPnt[3], double voxelSize, bool*** grids);
double CalcFractalDimension(const double* pntCoords, int nPnts,
    const int* triangle3Vers, int numTriangles, double porosity, int genus, bool bProgressBar)
{
    if (nPnts < 3)
        return -1.0;
    double boxSizes[] = { 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0};

    int numBoxSizes = sizeof(boxSizes) / sizeof(boxSizes[0]);
    int i, j, k;
    double boxSize;

    //step1: calculate model's bounding box
    double minPnt[3], maxPnt[3];
    minPnt[0] = minPnt[1] = minPnt[2] = std::numeric_limits<double>::max();
    maxPnt[0] = maxPnt[1] = maxPnt[2] = -std::numeric_limits<double>::max();
    for (i = 0, j = 0; i < nPnts; ++i, j += 3)
    {
        if (pntCoords[j] < minPnt[0])
            minPnt[0] = pntCoords[j];
        if (pntCoords[j + 1] < minPnt[1])
            minPnt[1] = pntCoords[j + 1];
        if (pntCoords[j + 2] < minPnt[2])
            minPnt[2] = pntCoords[j + 2];
        if (pntCoords[j] > maxPnt[0])
            maxPnt[0] = pntCoords[j];
        if (pntCoords[j + 1] > maxPnt[1])
            maxPnt[1] = pntCoords[j + 1];
        if (pntCoords[j + 2] > maxPnt[2])
            maxPnt[2] = pntCoords[j + 2];
    }
    double dx = maxPnt[0] - minPnt[0];
    double dy = maxPnt[1] - minPnt[1];
    double dz = maxPnt[2] - minPnt[2];
    double maxD = std::max(std::max(dx, dy), dz);
    if (maxD < 1.0e-7)
        return -1.0;

    //step2: uniformly scale raw model to {0,0,0} - {UNIFOMR_LEN, UNIFOMR_LEN, UNIFOMR_LEN};
    const double UNIFOMR_LEN = 100.0;
    std::vector<double> scaledPntCoords(pntCoords, pntCoords + 3 * nPnts);
    double scale = UNIFOMR_LEN / maxD;
    for (i = 0, j = 0; i < nPnts; ++i, j += 3)
    {
        scaledPntCoords[j] = (scaledPntCoords[j] - minPnt[0]) * scale;
        scaledPntCoords[j + 1] = (scaledPntCoords[j + 1] - minPnt[1]) * scale;
        scaledPntCoords[j + 2] = (scaledPntCoords[j + 2] - minPnt[2]) * scale;
    }
    minPnt[0] = minPnt[1] = minPnt[2] = 0;

    //step3: calculate taken grids number in different boxes sizes
    const double* coords = scaledPntCoords.data();
    long long boxCounts;
    std::vector<double> logBoxSizes, logBoxCounts;
    long long step = numBoxSizes * (numTriangles / 100), curStep = 0;
    bool*** grids;

    for (i = 0; i < numBoxSizes; ++i)
    {
        boxSize = boxSizes[i];

        int maxGrids = UNIFOMR_LEN / boxSize + 3;
        grids = new bool** [maxGrids];
        for (j = 0; j < maxGrids; ++j)
        {
            grids[j] = new bool* [maxGrids];
            for (k = 0; k < maxGrids; ++k)
            {
                grids[j][k] = new bool[maxGrids];
                memset(grids[j][k], 0, sizeof(bool) * maxGrids);
            }
        }
        double weight = AdjustWeight(porosity, genus);
        
        boxCounts = 0;
        long long totalBoxCounts = 0;

        for (j = 0; j < numTriangles; ++j)
            totalBoxCounts += FillTriangle(coords + 3 * triangle3Vers[3 * j], coords + 3 * triangle3Vers[3 * j + 1], coords + 3 * triangle3Vers[3 * j + 2],
                minPnt, boxSize, grids);
        boxCounts = totalBoxCounts * weight;

        for (j = 0; j < maxGrids; ++j)
        {
            for (k = 0; k < maxGrids; ++k)
                delete[] grids[j][k];
            delete[] grids[j];
        }
        delete[] grids;

        //exception
        if (boxCounts == 0)
            continue;

        logBoxSizes.push_back(log(boxSize));
        logBoxCounts.push_back(log(boxCounts + 0.0));
    }
    std::vector<std::array<double, 2>> graphPnts(logBoxSizes.size());
    for (i = 0; i < graphPnts.size(); ++i)
        graphPnts[i] = { logBoxSizes[i], logBoxCounts[i] };

    if (graphPnts.size() < 2)
        return -1;
    double a, b;
    bool bVertical;
    double aveError, maxError;

    if (!Geom_Fit_Line2D_LeastSquare(graphPnts, a, b, bVertical, aveError, maxError))
        return -1.0;
    return -a;
}

double mTranslate::AdjustWeight(double porosity, int genus)
{
    double porosityFactor = 0.0035;
    double genusFactor = 0.0031;
    if (genus < 0)
        genus = 0;

    double weight = 1 - porosity * porosityFactor + genus * genusFactor;

    if (weight < 0.0) {
        weight = 1.0;
    }

    return weight;
}

bool Geom_Fit_Line2D_LeastSquare(const std::vector<std::array<double, 2>>& points, double& a, double& b, bool& bVertical, double& aveError, double& maxError) {
    if (points.size() < 2)
        return false;

    double sum_x = 0.0, sum_y = 0.0;
    for (const auto& point : points) {
        sum_x += point[0];
        sum_y += point[1];
    }
    double mean_x = sum_x / points.size();
    double mean_y = sum_y / points.size();

    double numerator = 0.0, denominator = 0.0;
    for (const auto& point : points) {
        numerator += (point[0] - mean_x) * (point[1] - mean_y);
        denominator += std::pow(point[0] - mean_x, 2);
    }

    a = numerator / denominator;
    b = mean_y - a * mean_x;

    bVertical = (std::abs(a) < 1e-10);

    aveError = 0.0;
    maxError = 0.0;
    for (const auto& point : points) {
        double error = std::abs(point[1] - (a * point[0] + b));
        aveError += error;
        maxError = std::max(maxError, error);
    }
    aveError /= points.size();

    return true;
}

bool IsPointInsideTriangle(const double p[3], const double a[3], const double b[3], const double c[3])
{
    double u, v, w;
    u = ((b[1] - c[1]) * (p[0] - c[0]) + (c[0] - b[0]) * (p[1] - c[1])) / ((b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1]));
    v = ((c[1] - a[1]) * (p[0] - c[0]) + (a[0] - c[0]) * (p[1] - c[1])) / ((b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1]));
    w = 1.0f - u - v;

    return (u >= 0.0f) && (v >= 0.0f) && (w >= 0.0f);
}

//update method :fill triangle voxels   
long long mTranslate::FillTriangle(const double ver1[3], const double ver2[3], const double ver3[3],
    const double boxMinPnt[3], double voxelSize, bool*** grids)
{
    long long voxelCount = 0;

    int x0 = static_cast<int>((std::min(std::min(ver1[0], ver2[0]), ver3[0]) - boxMinPnt[0]) / voxelSize);
    int x1 = static_cast<int>((std::max(std::max(ver1[0], ver2[0]), ver3[0]) - boxMinPnt[0]) / voxelSize) + 1;
    int y0 = static_cast<int>((std::min(std::min(ver1[1], ver2[1]), ver3[1]) - boxMinPnt[1]) / voxelSize);
    int y1 = static_cast<int>((std::max(std::max(ver1[1], ver2[1]), ver3[1]) - boxMinPnt[1]) / voxelSize) + 1;
    int z0 = static_cast<int>((std::min(std::min(ver1[2], ver2[2]), ver3[2]) - boxMinPnt[2]) / voxelSize);
    int z1 = static_cast<int>((std::max(std::max(ver1[2], ver2[2]), ver3[2]) - boxMinPnt[2]) / voxelSize) + 1;

    std::mutex mtx;
#pragma omp parallel for collapse(3) reduction(+:voxelCount)
    for (int x = x0; x < x1; ++x)
    {
        for (int y = y0; y < y1; ++y)
        {
            for (int z = z0; z < z1; ++z)
            {
                double voxelPos[3] = { boxMinPnt[0] + (x + 0.5) * voxelSize,
                                      boxMinPnt[1] + (y + 0.5) * voxelSize,
                                      boxMinPnt[2] + (z + 0.5) * voxelSize };

                if (IsPointInsideTriangle(voxelPos, ver1, ver2, ver3))
                {
                    if (!grids[x][y][z])
                    {
#pragma omp critical
                        {
                            grids[x][y][z] = true;
                            voxelCount++;
                        }
                    }
                }
            }
        }
    }
    return voxelCount;
}

struct PairHash {
    template <typename T1, typename T2>
    std::size_t operator () (const std::pair<T1, T2>& pair) const {
        return std::hash<T1>{}(pair.first) ^ std::hash<T2>{}(pair.second);
    }
};

int mTranslate::CalculateGenus(int partID)
{
    std::vector<my::render::ViewData> viewdatas;
    model::AppImp::Instance()->_model->getPartById(partID)->getViewData(viewdatas);

    int numGenus = 0;

    for (const auto& viewdata : viewdatas) {
        std::vector<std::vector<double>> vertices;
        std::vector<std::vector<int>> faces;
        std::unordered_set<std::pair<int, int>, PairHash> edges; 

        if (!viewdata.isValid()) {
            Messager::Msg(Messager::eWarning, "Invalid part or no mesh data in this part!");
            return -1.0;
        }
        if (viewdata.wireframe)
            continue;

        int numVertices = viewdata.nodes.size() / 3;
        
        int numFaces = viewdata.trias.size() / 3;

        std::set<int> vers;

        for (size_t i = 0; i < viewdata.trias.size(); i += 3) {
            std::vector<int> face(viewdata.trias.begin() + i, viewdata.trias.begin() + i + 3);
            faces.push_back(face);
            for (size_t j = 0; j < face.size(); ++j) {
                int v1 = face[j];
                int v2 = face[(j + 1) % face.size()];
                vers.insert(v1);
                vers.insert(v2);
                std::pair<int, int> edge(std::min(v1, v2), std::max(v1, v2));

                if (edges.find(edge) == edges.end()) {
                    edges.insert(edge);
                }
            }
        }
        
        int numEdges = edges.size(); 
        
        numGenus = (2 - (numVertices - numEdges + numFaces)) / 2;
    }

    return numGenus;
}

double CalculateTriangleArea(const std::vector<double>& v1,
    const std::vector<double>& v2,
    const std::vector<double>& v3) {
    std::vector<double> edge1 = { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
    std::vector<double> edge2 = { v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2] };

    double cross_product_norm = sqrt(pow(edge1[1] * edge2[2] - edge1[2] * edge2[1], 2) +
        pow(edge1[2] * edge2[0] - edge1[0] * edge2[2], 2) +
        pow(edge1[0] * edge2[1] - edge1[1] * edge2[0], 2));

    double area = cross_product_norm / 2.0;

    return area;
}

static double CalculateVolume(const std::vector<std::vector<double>>& vertices, const std::vector<std::vector<int>>& triangles) {
    double volume = 0.0;

    for (const auto& triangle : triangles) {
        const std::vector<double>& v1 = vertices[triangle[0]];
        const std::vector<double>& v2 = vertices[triangle[1]];
        const std::vector<double>& v3 = vertices[triangle[2]];

        mTranslate obj;
        double tetrahedron_volume = obj.CalculateTriangleVolume(v1, v2, v3);

        volume += tetrahedron_volume;
    }

    volume = std::abs(volume);

    return volume;
}

double CalculateTriangleVolume (const std::vector<double>& v1,
    const std::vector<double>& v2,
    const std::vector<double>& v3) {
    double determinant = (v2[0] - v1[0]) * ((v3[1] - v1[1]) * v1[2] - (v3[2] - v1[2]) * v1[1]) -
        (v2[1] - v1[1]) * ((v3[0] - v1[0]) * v1[2] - (v3[2] - v1[2]) * v1[0]) +
        (v2[2] - v1[2]) * ((v3[0] - v1[0]) * v1[1] - (v3[1] - v1[1]) * v1[0]);

    double volume = determinant / 6.0;

    return volume;
}

double mTranslate::BoundingBoxVolume(const std::vector<std::vector<double>>& vertices) {
    if (vertices.empty()) {
        return 0.0;
    }
    double minX = vertices[0][0];
    double minY = vertices[0][1];
    double minZ = vertices[0][2];
    double maxX = vertices[0][0];
    double maxY = vertices[0][1];
    double maxZ = vertices[0][2];

    for (const auto& vertex : vertices) {
        minX = std::min(minX, vertex[0]);
        minY = std::min(minY, vertex[1]);
        minZ = std::min(minZ, vertex[2]);
        maxX = std::max(maxX, vertex[0]);
        maxY = std::max(maxY, vertex[1]);
        maxZ = std::max(maxZ, vertex[2]);
    }

    double volume = (maxX - minX) * (maxY - minY) * (maxZ - minZ);
    return volume;
}

double CalculateVolumeRatio(int partID) {
    std::vector<my::render::ViewData> viewdatas;
    model::AppImp::Instance()->_model->getPartById(partID)->getViewData(viewdatas);

    double totalVolume = 0.0;
    std::vector<std::vector<double>> allVertices;
    mTranslate obj;

    for (const auto& viewdata : viewdatas) {
        if (!viewdata.isValid()) {
            Messager::Msg(Messager::eWarning, "Invalid part or no mesh data in this part!");
            return -1.0;
        }

        for (size_t i = 0; i < viewdata.trias.size(); i += 3) {
            std::vector<double> v1(viewdata.nodes.begin() + 3 * viewdata.trias[i],
                viewdata.nodes.begin() + 3 * viewdata.trias[i] + 3);
            std::vector<double> v2(viewdata.nodes.begin() + 3 * viewdata.trias[i + 1],
                viewdata.nodes.begin() + 3 * viewdata.trias[i + 1] + 3);
            std::vector<double> v3(viewdata.nodes.begin() + 3 * viewdata.trias[i + 2],
                viewdata.nodes.begin() + 3 * viewdata.trias[i + 2] + 3);

            double volumeTriangle = obj.CalculateTriangleVolume(v1, v2, v3);
            totalVolume += volumeTriangle;

            allVertices.push_back(v1);
            allVertices.push_back(v2);
            allVertices.push_back(v3);
        }
    }
    double totalBoxVolume = obj.BoundingBoxVolume(allVertices);

    double volumeRatio = totalVolume / totalBoxVolume;

    double minVal = 0.0;
    double maxVal = 1.0;
    double percentageOfVR = (volumeRatio - minVal) / (maxVal - minVal) * 100.0;

    return percentageOfVR;
}


