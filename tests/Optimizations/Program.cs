using System.Reflection;
using System.Diagnostics;
using System.Globalization;
using System.Collections.Immutable;
using System.Security.Cryptography;
using System.Text.Json;

class Harness
{
    static Type p = null!;
    static object? Get(string n) => p.GetField(n)!.GetValue(null);
    static void Set(string n, object v) => p.GetField(n)!.SetValue(null,v);
    static object? Call(string n, params object[] args) => p.GetMethod(n,BindingFlags.Static|BindingFlags.Public|BindingFlags.NonPublic)!.Invoke(null,args);
    static void Main(string[] args)
    {
        CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
        p=Assembly.LoadFrom(Path.GetFullPath(args[0])).GetType("GRAMM_2001.Program")!;
        int nx=80,ny=70,nz=32;
        Set("NX",nx);Set("NY",ny);Set("NZ",nz);Set("pOptions",new ParallelOptions{MaxDegreeOfParallelism=1});
        GC.Collect(); long before=GC.GetTotalMemory(true); var sw=Stopwatch.StartNew(); Call("Define_Arrays"); sw.Stop();
        long memory=GC.GetTotalMemory(true)-before;
        foreach(string name in new[]{"RHO","AP0"})
        {
            var a=(float[][][])Get(name)!; for(int i=0;i<a.Length;i++) for(int j=0;j<a[i].Length;j++) for(int k=0;k<a[i][j].Length;k++) a[i][j][k]=1f+(i+j+k)%17*.01f;
        }
        foreach(string name in new[]{"U1N","U2N","V1N","V2N","W1N","W2N","DP","DPX","DPY","DPZ"})
        {
            var a=(double[][][])Get(name)!; int seed=name.Sum(c=>(int)c);
            for(int i=0;i<a.Length;i++) for(int j=0;j<a[i].Length;j++) for(int k=0;k<a[i][j].Length;k++) a[i][j][k]=Math.Sin(seed+i*3+j*7+k)*.01;
        }
        foreach(string name in new[]{"AREAImm","AREAXImm","AREAYImm","AREAZImm","AREAZXImm","AREAZYImm","AREAXYZImm"})
        {
            var a=(ImmutableArray<float>[][])Get(name)!;
            for(int i=0;i<a.Length;i++) for(int j=0;j<a[i].Length;j++) a[i][j]=Enumerable.Range(0,nz+2).Select(k=>1f+(i+j+k)%9*.03f).ToImmutableArray();
        }
        Set("ICPN",false);
        var timing=new List<double>();
        for(int repeat=0;repeat<8;repeat++) {sw.Restart();Call("CALCPR_calculate",nx,ny,nz);sw.Stop();timing.Add(sw.Elapsed.TotalMilliseconds);}
        var hashes=new Dictionary<string,string>();
        foreach(string name in new[]{"U1N","U2N","V1N","V2N","W1N","W2N","UN","VN","WN","SUX","SUY","SUZ","SUXYZ","DP","DPX","DPY","DPZ"})
        {
            using var hash=IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
            foreach(Array row in (Array)Get(name)!) foreach(Array col in row) {var b=new byte[Buffer.ByteLength(col)];Buffer.BlockCopy(col,0,b,0,b.Length);hash.AppendData(b);}
            hashes[name]=Convert.ToHexString(hash.GetHashAndReset());
        }
        Call("ClearArrays");
        File.WriteAllText(args[1],JsonSerializer.Serialize(new{status="pass",nx,ny,nz,managed_bytes=memory,median_ms=timing.Skip(1).Order().ElementAt(3),samples_ms=timing,hashes,SUMG=Get("SUMG")},new JsonSerializerOptions{WriteIndented=true}));
        Console.WriteLine("PASS GRAMM arrays");
    }
}
