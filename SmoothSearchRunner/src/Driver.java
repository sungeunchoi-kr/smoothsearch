import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.nio.charset.StandardCharsets;
import java.util.Scanner;

public class Driver {

	private static String baseUrl = "http://162.243.237.110:8718";
	private static URL getFreeBlockAddressUrl;
	private static URL postSolvedBlockUrlBase;
	
	private static void initialize() {
		try {
			getFreeBlockAddressUrl = new URL(baseUrl + "/v1/free-blocks");
			postSolvedBlockUrlBase = new URL(baseUrl + "/v1/searched-blocks/");
		} catch (MalformedURLException e) {
			e.printStackTrace();
		}
	}
	
	public static String getBlockAddress() {
		try {
			final URLConnection conn = getFreeBlockAddressUrl.openConnection();
			final InputStream response = conn.getInputStream();		
			final Scanner scanner = new Scanner(response);
		    final String blockAddress = scanner.useDelimiter("\\A").next();
		    scanner.close();
		    response.close();
		    String t = blockAddress.substring(1, blockAddress.length()-1);
			return t;
		} catch (IOException e) {
			e.printStackTrace();
			System.out.println("EXCEPTION: " + e.toString());
			return null;
		}
	}
	
	public static boolean postSolvedBlock(String blockAddress, String payload) {
		try {
			final URLConnection connection = new URL(postSolvedBlockUrlBase, blockAddress).openConnection();
			connection.setDoOutput(true);
			//connection.setRequestProperty("Accept-Charset", charset);
			//connection.setRequestProperty("Content-Type", "application/x-www-form-urlencoded;charset=" + charset);
	
			final OutputStream output = connection.getOutputStream();
			output.write(payload.getBytes(StandardCharsets.UTF_8));
			final InputStream response = connection.getInputStream();
			final Scanner scanner = new Scanner(response);
		    final String serverResponse = scanner.useDelimiter("\\A").next();
		    System.out.println("serverResponse=" + serverResponse);
		    scanner.close();
			output.close();
			return true;
		} catch (IOException e) {
			e.printStackTrace();
			System.out.println("EXCEPTION: " + e.toString());
			return false;
		}
	}
	
	public static void main(String[] args) {
		System.out.println("program started.");
		initialize();
		try {
			while (true) {
				System.out.println("requesting free block address...");
				final String addr = getBlockAddress();
				System.out.println("got block address " + addr + ".");
				Process p = new ProcessBuilder("./run.exe", addr).start();
				InputStream procResult = p.getInputStream();
				final Scanner scanner = new Scanner(procResult);
			    final String procResultString = scanner.useDelimiter("\\A").next();
			    scanner.close();
			    procResult.close();
			    System.out.println("procResultString: " + procResultString);
				postSolvedBlock(addr, procResultString);
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.out.println("EXCEPTION: " + e.toString());
		}
	}
}
