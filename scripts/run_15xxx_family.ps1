param(
    [int]$MaxParallel = 4,
    [int]$PollSeconds = 5,
    [int]$MaxRetries = 3,
    [int]$MatlabGraceSeconds = 180,   # tempo concesso a MATLAB per avviarsi prima di controllarlo
    [int]$StallMinutes = 30           # log FAST.Farm fermo per piu' di cosi' --> caso bloccato
)

$FastfarmDirWsl = "/mnt/d/zero_model/fastfarm_simulation"
$PyDirWsl       = "/mnt/d/zero_model/controller_floris"
$MatlabDirWin   = "D:\zero_model\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3"
$CondaEnv       = "rosco-env"
$MatlabExe      = "C:\Program Files\MATLAB\R2024a\bin\matlab.exe"
$LogDir         = "D:\zero_model\logs_15xxx"
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null

$BaselinePorts = @(15006, 15012, 15106, 15112, 15206, 15212)

# Porte da saltare (casi già completati). Vuota: nuova campagna con TMax=4200, DT_High=0.25,
# inflow seed_wind8_DELNN_new_seed --> si rilanciano tutte le 36 porte.
$AlreadyDone = @()
$AllPorts = ((15001..15012) + (15101..15112) + (15201..15212)) | Where-Object { $_ -notin $AlreadyDone }

$queue = [System.Collections.Generic.Queue[int]]::new()
$AllPorts | ForEach-Object { $queue.Enqueue($_) }
$active = @{}          # port -> pscustomobject{FF; Py; Matlab; IsBaseline}
$retryCount = @{}       # port -> numero di tentativi falliti finora
$failedForGood = @()    # porte che hanno esaurito i tentativi

function Start-Case([int]$port) {
    $isBaseline = $BaselinePorts -contains $port
    $logFF = Join-Path $LogDir "$port`_fastfarm.log"

    # Il controller Python serve SEMPRE (anche i baseline hanno YCMode=5 e vogliono un client ZMQ
    # connesso), ma per i baseline è la variante "greedy" senza MATLAB -- vedi FLORIDyn_<port>.py.
    $logPy = Join-Path $LogDir "$port`_python.log"
    $pyBin = "/home/tente/miniforge3/envs/$CondaEnv/bin/python"
    $pyCmd = "cd $PyDirWsl && $pyBin FLORIDyn_$port.py"
    $pyArgStr = "bash -ic `"$pyCmd`""
    $py = $null
    try {
        $py = Start-Process -FilePath "wsl.exe" -ArgumentList $pyArgStr `
              -RedirectStandardOutput $logPy -RedirectStandardError "$logPy.err" `
              -WindowStyle Hidden -PassThru -ErrorAction Stop
    } catch {
        Write-Host "[$(Get-Date -Format 'HH:mm:ss')] $port  Start-Process del controller Python fallito: $_"
    }

    $ml = $null
    if (-not $isBaseline) {
        $logMl = Join-Path $LogDir "$port`_matlab.log"
        # try/catch: se lo script va in errore MATLAB deve uscire (senza, resta aperto e inattivo e
        # FAST.Farm aspetta per sempre). -logfile perche' su Windows MATLAB non scrive su stdout.
        $mlExpr = "try, run('$MatlabDirWin\FLORIDyn_main_ZMQ_$port.m'); catch e, disp(getReport(e)); end; exit;"
        $mlArgStr = "-nosplash -nodesktop -minimize -logfile `"$logMl`" -r `"$mlExpr`""
        $ml = Start-Process -FilePath $MatlabExe -ArgumentList $mlArgStr `
              -WindowStyle Hidden -PassThru
    }

    Start-Sleep -Seconds 3   # dai tempo a controller (+matlab) prima di lanciare FAST.Farm

    $ffCmd = "~/bin/FF42 $FastfarmDirWsl/$port.fstf"
    $ffArgStr = "bash -ic `"$ffCmd`""
    $ff = Start-Process -FilePath "wsl.exe" -ArgumentList $ffArgStr `
          -RedirectStandardOutput $logFF -RedirectStandardError "$logFF.err" `
          -WindowStyle Hidden -PassThru

    Write-Host "[$(Get-Date -Format 'HH:mm:ss')] START $port (baseline=$isBaseline)"
    $now = Get-Date
    return [pscustomobject]@{ FF = $ff; Py = $py; Matlab = $ml; IsBaseline = $isBaseline;
                              Started = $now; LogFF = $logFF; LogLen = -1; LastGrowth = $now }
}

# Il PID di Start-Process per matlab.exe e' quello del launcher, non della sessione reale:
# cerco la sessione per riga di comando (stesso criterio usato in Stop-Case).
function Test-MatlabAlive([int]$port) {
    $p = Get-CimInstance Win32_Process -Filter "Name='MATLAB.exe'" -ErrorAction SilentlyContinue |
         Where-Object { $_.CommandLine -like "*FLORIDyn_main_ZMQ_$port.m*" }
    return [bool]$p
}

# Lunghezza del log letta da un handle aperto: i metadati di cartella di un file ancora
# aperto in scrittura possono non aggiornarsi.
function Get-LogLength([string]$path) {
    try {
        $fs = [System.IO.File]::Open($path, 'Open', 'Read', 'ReadWrite')
        $len = $fs.Length
        $fs.Close()
        return $len
    } catch { return -1 }
}

function Stop-Case([int]$port, [pscustomobject]$case) {
    foreach ($p in @($case.FF, $case.Py, $case.Matlab)) {
        if ($p -and -not $p.HasExited) {
            try { Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue } catch {}
        }
    }
    # Fallback robusto: il PID catturato da Start-Process per matlab.exe a volte non e' quello
    # della sessione reale (si e' visto restare vive per giorni sessioni "chiuse" con successo).
    # Cerca e chiudi per riga di comando invece di fidarti solo dell'handle.
    Get-CimInstance Win32_Process -Filter "Name='MATLAB.exe'" -ErrorAction SilentlyContinue |
        Where-Object { $_.CommandLine -like "*FLORIDyn_main_ZMQ_$port.m*" } |
        ForEach-Object { try { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue } catch {} }
}

function Retry-OrGiveUp([int]$port) {
    $n = 1
    if ($retryCount.ContainsKey($port)) { $n = $retryCount[$port] + 1 }
    $retryCount[$port] = $n
    if ($n -le $MaxRetries) {
        Write-Host "[$(Get-Date -Format 'HH:mm:ss')] RETRY $port (tentativo $n/$MaxRetries)"
        $queue.Enqueue($port)
    } else {
        Write-Host "[$(Get-Date -Format 'HH:mm:ss')] FALLITO DEFINITIVAMENTE $port dopo $MaxRetries tentativi"
        $script:failedForGood += $port
    }
}

# --- loop principale: riempi gli slot fino a MaxParallel, poi polla ---
while ($active.Count -gt 0 -or $queue.Count -gt 0) {
    while ($active.Count -lt $MaxParallel -and $queue.Count -gt 0) {
        $port = $queue.Dequeue()
        $active[$port] = Start-Case $port
    }

    Start-Sleep -Seconds $PollSeconds

    foreach ($port in @($active.Keys)) {
        $case = $active[$port]

        if ($case.FF.HasExited) {
            $code = $case.FF.ExitCode
            Stop-Case $port $case
            Write-Host "[$(Get-Date -Format 'HH:mm:ss')] END   $port  FAST.Farm exit=$code"
            $active.Remove($port)
            continue
        }

        # Fallimento precoce: il controller Python e' morto (o non e' mai partito, es. Start-Process
        # fallito) ma FAST.Farm e' ancora vivo --> e' orfano e restera' bloccato per sempre in
        # attesa di ZMQ. Meglio abortire subito e ritentare, invece di scoprirlo dopo giorni.
        if ((-not $case.Py) -or $case.Py.HasExited) {
            Write-Host "[$(Get-Date -Format 'HH:mm:ss')] ABORT $port  (controller Python morto/mai partito, FAST.Farm orfano -- termino e rimetto in coda)"
            Stop-Case $port $case
            $active.Remove($port)
            Retry-OrGiveUp $port
            continue
        }

        # MATLAB morto (crash o errore nello script): FAST.Farm resterebbe bloccato in attesa
        # della risposta di FLORIDyn. Controllo solo dopo il tempo di avvio di MATLAB.
        if ((-not $case.IsBaseline) -and
            ((Get-Date) - $case.Started).TotalSeconds -gt $MatlabGraceSeconds -and
            (-not (Test-MatlabAlive $port))) {
            Write-Host "[$(Get-Date -Format 'HH:mm:ss')] ABORT $port  (MATLAB morto -- vedi $port`_matlab.log -- termino e rimetto in coda)"
            Stop-Case $port $case
            $active.Remove($port)
            Retry-OrGiveUp $port
            continue
        }

        # Stallo generico: il log di FAST.Farm non cresce da StallMinutes (qualunque processo
        # si sia bloccato, la simulazione non avanza).
        $len = Get-LogLength $case.LogFF
        if ($len -ne $case.LogLen) {
            $case.LogLen = $len
            $case.LastGrowth = Get-Date
        } elseif (((Get-Date) - $case.LastGrowth).TotalMinutes -gt $StallMinutes) {
            Write-Host "[$(Get-Date -Format 'HH:mm:ss')] ABORT $port  (FAST.Farm fermo da $StallMinutes min -- termino e rimetto in coda)"
            Stop-Case $port $case
            $active.Remove($port)
            Retry-OrGiveUp $port
            continue
        }
    }
}

Write-Host "Coda esaurita. Log in $LogDir"
if ($failedForGood.Count -gt 0) {
    Write-Host "Porte fallite definitivamente dopo $MaxRetries tentativi: $($failedForGood -join ', ')"
} else {
    Write-Host "Tutte le porte completate con successo."
}
