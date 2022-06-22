import os
import subprocess
subprocess.run('pipenv install', shell=True)
import telegram
import json
import waitress
import main

BOT_TOKEN = '5474046498:AAFYBH-DG_Ob4sDatPAyb9ap6Mev4s9aQNU'
BOT_URL = f"https://api.telegram.org/bot{BOT_TOKEN}"  
TMP_JSON = "./options.json"

# download ngrok
subprocess.run('curl -s https://ngrok-agent.s3.amazonaws.com/ngrok.asc | tee \
              /etc/apt/trusted.gpg.d/ngrok.asc >/dev/null && echo "deb \
              https://ngrok-agent.s3.amazonaws.com/ buster main" | tee \
              /etc/apt/sources.list.d/ngrok.list && apt update && apt \
              install ngrok', shell=True, check=True)


def receive_server_options():
    os.system(f'curl localhost:4040/api/tunnels > {TMP_JSON}')
    try: 
        with open(TMP_JSON) as f:
            options = json.load(f)
        os.remove(TMP_JSON)
    except BaseException:
        return {'tunnels': []}
    return options
      

bot = telegram.Bot(token=BOT_TOKEN)

if __name__ == '__main__':
    cmd = ['ngrok', 'http', '8443']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
 
    options = receive_server_options()
    while not options['tunnels']:
         options = receive_server_options()

    ngrok_url = options['tunnels'][0]['public_url']

    bot.setWebhook('{URL}/{HOOK}'.format(URL=ngrok_url, HOOK=BOT_TOKEN))
    print(ngrok_url)
    print(bot)
    
    waitress.serve(main.app,
                   host="0.0.0.0",
                   port=8443)
     
    print(process.stdout.read())
    process.stdout.close()