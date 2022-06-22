import telegram
from flask import Flask, request
from configure_and_run import bot, BOT_TOKEN, BOT_URL 

app = Flask(__name__)

@app.route(f'/{BOT_TOKEN}', methods=['POST'])
def webhook_post():

    received_message = request.json['message']
    chat_id = received_message['chat']['id']
    print(received_message)
    if 'text' in received_message:
        received_text = received_message['text']
        bot.sendMessage(chat_id, received_text)
    elif 'photo' in received_message:
        bot.sendPhoto(chat_id, received_message['photo'][0]['file_id'])  

    elif ('document' in received_message and 
          received_message['document']['file_name'].endswith(('.png', '.jpg'))):
        bot.sendDocument(chat_id, received_message['document']['file_id'])    
    else:
        bot.sendMessage(chat_id, 'Bot can answer only on text messages, photos and files .png .jpg!')
    return 'Ok'        

