import java.util.ArrayList;
import java.util.Collections;

public class Main {
    public static void printArraylistIntegerArray(ArrayList<Integer[]> old){
        String newStr = "";
        for (int i = 0; i < old.size(); i++){
            Integer[] array = old.get(i);
            for (int j = 0; j < array.length; j++){
                if (j == 0){
                    newStr += "[";
                }
                newStr += array[j] + ", ";
                if (j == array.length-1){
                    newStr = newStr.substring(0, newStr.length()-2);
                    newStr += "], ";
                }
            }
            if (i == old.size()-1){
                newStr = newStr.substring(0, newStr.length()-2);
            }
        }
        System.out.println(newStr);
    }

    public static void printTwoDimentionalArray(Integer[][] old){
        String newStr = "";
        for (int i = 0; i < old.length; i++){
            Integer[] array = old[i];
            for (int j = 0; j < array.length; j++){
                if (j == 0){
                    newStr += "[";
                }
                newStr += array[j] + ", ";
                if (j == array.length-1){
                    newStr = newStr.substring(0, newStr.length()-2);
                    newStr += "], ";
                }
            }
            if (i == old.length-1){
                newStr = newStr.substring(0, newStr.length()-2);
            }
        }
        System.out.println(newStr);
    }

    public static int swapCost(int firstPointXCoordinate, int firstPointYCoordinate, int secondPointXCoordinate, int secondPointYCoordinate,
                               int thirdPointXCoordinate, int thirdPointYCoordinate, int fourthPointXCoordinate, int fourthPointYCoordinate){
        // swapCost( (firstPoint, secondPoint), (thirdPoint, fourthPoint) ) = || firstPoint - fourthPoint || + || thirdPoint- secondPoint || - || firstPoint - secondPoint || - || thirdPoint - fourthPoint ||

        // || firstPoint - fourthPoint ||
        int firstPoint_fourthPoint = (int) Math.abs(euclideanDistance(firstPointXCoordinate, firstPointYCoordinate, fourthPointXCoordinate, fourthPointYCoordinate));

        // || thirdPoint- secondPoint ||
        int thirdPoint_secondPoint = (int) Math.abs(euclideanDistance(thirdPointXCoordinate, thirdPointYCoordinate, secondPointXCoordinate, secondPointYCoordinate));

        // || firstPoint - secondPoint ||
        int firstPoint_secondPoint = (int) Math.abs(euclideanDistance(firstPointXCoordinate, firstPointYCoordinate, secondPointXCoordinate, secondPointYCoordinate));

        // || thirdPoint - fourthPoint ||
        int thirdPoint_fourthPoint = (int) Math.abs(euclideanDistance(thirdPointXCoordinate, thirdPointYCoordinate, fourthPointXCoordinate, fourthPointYCoordinate));

        return firstPoint_fourthPoint + thirdPoint_secondPoint - firstPoint_secondPoint - thirdPoint_fourthPoint;
    }

    public static double euclideanDistance(double x1, double y1, double x2, double y2){
        return Math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }

    public static void main(String[] args) {
        ArrayList<Integer[]> globalTpsPathMatrix = new ArrayList<>();
        ArrayList<Integer[]> blockTpsPathMatrix = new ArrayList<>();

        Integer[] coordinate1 = {1,2};
        Integer[] coordinate2 = {35,4};
        Integer[] coordinate3 = {7,55};
        Integer[] coordinate4 = {44,4};
        Integer[] coordinate5 = {14,23};
        Integer[] coordinate6 = {16,2};
        Integer[] coordinate7 = {14,22};
        Integer[] coordinate8 = {16,27};
        Integer[] coordinate9 = {12,25};
        Integer[] coordinate10 = {1,26};
        globalTpsPathMatrix.add(coordinate1);
        globalTpsPathMatrix.add(coordinate2);
        globalTpsPathMatrix.add(coordinate3);
        globalTpsPathMatrix.add(coordinate4);
        globalTpsPathMatrix.add(coordinate5);
        globalTpsPathMatrix.add(coordinate6);
        globalTpsPathMatrix.add(coordinate7);
        globalTpsPathMatrix.add(coordinate8);
        globalTpsPathMatrix.add(coordinate9);
        globalTpsPathMatrix.add(coordinate10);

        Integer[] coordinate11 = {13,25};
        Integer[] coordinate12 = {13,22};
        Integer[] coordinate13 = {41,27};
        Integer[] coordinate14 = {13,24};
        Integer[] coordinate15 = {15,28};
        blockTpsPathMatrix.add(coordinate11);
        blockTpsPathMatrix.add(coordinate12);
        blockTpsPathMatrix.add(coordinate13);
        blockTpsPathMatrix.add(coordinate14);
        blockTpsPathMatrix.add(coordinate15);

        System.out.print("Global Tsp Matrix Coordinates: \t");
        printArraylistIntegerArray(globalTpsPathMatrix);
        System.out.print("Block Tsp Matrix Coordinates: \t");
        printArraylistIntegerArray(blockTpsPathMatrix);
        System.out.print("Stitching Coordinates: \t\t\t");
        printTwoDimentionalArray(findCoordinates(globalTpsPathMatrix, blockTpsPathMatrix)); // first two blockCoordinates, last two GlobalCoordinates
        System.out.print("After Stitching: \t\t\t\t");
        printArraylistIntegerArray(stitchingCoordinates(globalTpsPathMatrix, blockTpsPathMatrix));


    }

    public static Integer[][] findCoordinates(ArrayList<Integer[]> globalTpsPathMatrix, ArrayList<Integer[]> blockTpsPathMatrix){
        int minimumSwapCost = (int) Double.POSITIVE_INFINITY;
        Integer[][] minimumSwapCostCoordinates = new Integer[4][2];

        for (int i = 0; i < blockTpsPathMatrix.size() - 1; i ++){

            Integer[] firstPoint = blockTpsPathMatrix.get(i);
            Integer[] secondPoint = blockTpsPathMatrix.get(i+1);

            int firstPointXCoordinate = blockTpsPathMatrix.get(i)[0];
            int firstPointYCoordinate = blockTpsPathMatrix.get(i)[1];

            int secondPointXCoordinate = blockTpsPathMatrix.get(i+1)[0];
            int secondPointYCoordinate = blockTpsPathMatrix.get(i+1)[1];

//            System.out.println("First Point: [" + String.valueOf(firstPointXCoordinate) + ", " + String.valueOf(firstPointYCoordinate) + "], Second Point: [" + String.valueOf(secondPointXCoordinate) + ", " + String.valueOf(secondPointYCoordinate) + "]" );

            for (int j = 0; j < globalTpsPathMatrix.size() -1; j ++){

                Integer[] thirdPoint = globalTpsPathMatrix.get(j);
                Integer[] fourthPoint = globalTpsPathMatrix.get(j+1);

                int thirdPointXCoordinate = globalTpsPathMatrix.get(j)[0];
                int thirdPointYCoordinate = globalTpsPathMatrix.get(j)[1];

                int fourthPointXCoordinate = globalTpsPathMatrix.get(j+1)[0];
                int fourthPointYCoordinate = globalTpsPathMatrix.get(j+1)[1];

                int swapCost = swapCost(firstPointXCoordinate, firstPointYCoordinate, secondPointXCoordinate, secondPointYCoordinate,
                        thirdPointXCoordinate, thirdPointYCoordinate, fourthPointXCoordinate, fourthPointYCoordinate);

                if (swapCost < minimumSwapCost){
                    minimumSwapCost = swapCost;
                    minimumSwapCostCoordinates[0] = firstPoint;
                    minimumSwapCostCoordinates[1] = secondPoint;
                    minimumSwapCostCoordinates[2] = thirdPoint;
                    minimumSwapCostCoordinates[3] = fourthPoint;
                }
            }
        }

        return minimumSwapCostCoordinates;
    }

    public static ArrayList<Integer[]> stitchingCoordinates(ArrayList<Integer[]> globalTpsPathMatrix, ArrayList<Integer[]> blockTpsPathMatrix){
        ArrayList<Integer[]> afterStitching = new ArrayList<>();

        Integer[][] coordinates = findCoordinates(globalTpsPathMatrix, blockTpsPathMatrix);

        // Block Coorodinates
        Integer[] firstPoint = coordinates[0];
        Integer[] secondPoint = coordinates[1];

        // Global Coorodinates
        Integer[] thirdPoint = coordinates[2];
        Integer[] fourthPoint = coordinates[3];

        // Block Index
        int firstPointIndex = blockTpsPathMatrix.indexOf(firstPoint);
        int secondPointIndex = blockTpsPathMatrix.indexOf(secondPoint);

        // Global Index
        int thirdPointIndex = globalTpsPathMatrix.indexOf(thirdPoint);
        int fourthPointIndex = globalTpsPathMatrix.indexOf(fourthPoint);

        if (globalTpsPathMatrix.size() == 0){
            afterStitching.addAll(blockTpsPathMatrix);
            return afterStitching;
        }

        else if (blockTpsPathMatrix.size() == 0){
            afterStitching.addAll(globalTpsPathMatrix);
            return afterStitching;
        }


        // First Leg -> before the global coordinates
        for (int i = 0; i <= thirdPointIndex; i++ ){
            afterStitching.add(globalTpsPathMatrix.get(i));
        }

        // Second Leg -> before the block coordinates
        ArrayList<Integer[]> temp = new ArrayList<>();
        for (int i = 0; i <= firstPointIndex; i++ ){
            temp.add(blockTpsPathMatrix.get(i));
        }
        Collections.reverse(temp);
        afterStitching.addAll(temp);

        // Third Leg -> After the block coordinates
        temp.clear();
        for (int i = secondPointIndex; i < blockTpsPathMatrix.size(); i++ ){
            temp.add(blockTpsPathMatrix.get(i));
        }
        Collections.reverse(temp);
        afterStitching.addAll(temp);

        // Third Leg -> After the global coordinates
        temp.clear();
        for (int i = fourthPointIndex; i < globalTpsPathMatrix.size(); i++ ){
            temp.add(globalTpsPathMatrix.get(i));
        }
        Collections.reverse(temp);
        afterStitching.addAll(temp);


        return afterStitching;
    }
}
