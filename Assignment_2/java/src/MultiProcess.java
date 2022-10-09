public class MultiProcess {
    public static void main(String[] args) {
        ProcessBuilder pb = new ProcessBuilder();
        pb.inheritIO(); // <-- passes IO from forked process.
        try {
            Process p = pb.start(); // <-- forkAndExec on Unix
            p.waitFor(); // <-- waits for the forked process to complete.
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
