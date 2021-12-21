package org.vmk.dep508.stream.iris;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;


public class App {
    public static void main(String[] args) throws IOException {
        App a = new App();
        a.test();
    }

    public void test() throws IOException {
        //load data from file iris.data
        List<Iris> irisList = Files.lines(Paths.get("iris.data"))
                .map(Iris::parse)
                .collect(Collectors.toList());
        IrisDataSetHelper helper = new IrisDataSetHelper(irisList);

        //get average sepal width
        Double avgSetalLength = helper.getAverage(Iris::getSepalWidth);
        System.out.println(avgSetalLength);

        //get average petal square - petal width multiplied on petal length
        Double avgPetalLength = helper.getAverage((x) -> x.getPetalLength() * x.getPetalWidth());
        System.out.println(avgPetalLength);

        //get average petal square for flowers with sepal width > 4
        Double avgPetalSquare = helper.getAverageWithFilter(x -> x.getSepalWidth() > 4,
                                                            x -> x.getPetalLength() * x.getPetalWidth());
        System.out.println(avgPetalSquare);

        //get flowers grouped by Petal size (Petal.SMALL, etc.)
        Map groupsByPetalSize = helper.groupBy(Iris::classifyByPatel);
        System.out.println(groupsByPetalSize);

        //get max sepal width for flowers grouped by species
        Map maxSepalWidthForGroupsBySpecies = helper.maxFromGroupedBy(Iris::getSpecies, Iris::getSepalWidth);
        System.out.println(maxSepalWidthForGroupsBySpecies);
    }

}

